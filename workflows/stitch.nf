#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for preparing lcWGS data for haplotype reconstruction
// using STITCH
// General flow:
// 1) Create channel of paired end or single end reads
// 1a) Concatenate reads from the same sample or strain if desired
// 2) Calculate quality statistics for sequencing
// 3) Generate read groups
// 4) BWA-mem alignment
// 5) Sort SAM files
// 6) Mark duplicates and gather alignment summary information
// 7) Collect .bam filenames in its own list
// 8) Generate other required input files for STITCH
// 9) Run STITCH
// 10) Aggregate summary statistics

// import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/stitch.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {TRIMMOMATIC_PE} from "${projectDir}/modules/utility_modules/trimmomatic"
include {FASTQC} from "${projectDir}/modules/utility_modules/fastqc"
include {MULTIQC} from "${projectDir}/modules/utility_modules/multiqc"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"
include {AGGREGATE_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_wgs"
include {CREATE_BAMLIST} from "${projectDir}/modules/utility_modules/create_bamlist"
include {CREATE_POSFILE} from "${projectDir}/modules/bcftools/create_posfile"
include {CREATE_POSFILE_DO} from "${projectDir}/modules/bcftools/create_posfile_DO"
include {RUN_STITCH} from "${projectDir}/modules/stitch/run_stitch"
include {RUN_STITCH_DO} from "${projectDir}/modules/stitch/run_stitch_DO"
include {STITCH_VCF_TO_TXT} from "${projectDir}/modules/stitch/vcf_to_sample_genos"
include {STITCH_TO_QTL} from "${projectDir}/modules/stitch/stitch_to_qtl2files"
include {STATS_MARKDOWN} from "${projectDir}/modules/utility_modules/render_stats_markdown"

// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// prepare reads channel
if (params.concat_lanes){
  if (params.read_type == 'PE'){
    read_ch = Channel
            .fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true, flat:true )
            .map { file, file1, file2 -> tuple(getLibraryId(file), file1, file2) }
            .groupTuple()
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}", checkExists:true, size:1 )
                .map { file, file1 -> tuple(getLibraryId(file), file1) }
                .groupTuple()
                .map{t-> [t[0], t[1].flatten()]}
  }
} else {
  if (params.read_type == 'PE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/${params.pattern}${params.extension}",checkExists:true )
  }
  else if (params.read_type == 'SE'){
    read_ch = Channel.fromFilePairs("${params.sample_folder}/*${params.extension}",checkExists:true, size:1 )
  }
}

// if channel is empty give error message and exit
read_ch.ifEmpty{ exit 1, "ERROR: No Files Found in Path: ${params.sample_folder} Matching Pattern: ${params.pattern}"}

chrs = Channel.of(1..19,'X','Y')

// main workflow
workflow STITCH {

  // Step 0: Concatenate Fastq files if required. 
  if (params.concat_lanes){
    if (params.read_type == 'PE'){
        CONCATENATE_READS_PE(read_ch)
        read_ch = CONCATENATE_READS_PE.out.concat_fastq
    } else if (params.read_type == 'SE'){
        CONCATENATE_READS_SE(read_ch)
        read_ch = CONCATENATE_READS_SE.out.concat_fastq
    }
  }

  // Run trimmomatic
  TRIMMOMATIC_PE(read_ch)

  // Calculate quality statistics for sequencing
  QUALITY_STATISTICS(TRIMMOMATIC_PE.out.trimmomatic)
  
  // Run fastqc on adapter trimmed reads
  FASTQC(TRIMMOMATIC_PE.out.to_fastqc)

  // Run multiqc
  fastqc_reports = FASTQC.out.to_multiqc.flatten().collect()
  MULTIQC(fastqc_reports)

  // Generate read groups
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "gatk")
  bwa_mem_mapping = QUALITY_STATISTICS.out.trimmed_fastq
                    .join(READ_GROUPS.out.read_groups)

  // BWA-mem alignment
  BWA_MEM(bwa_mem_mapping)

  // Sort SAM files
  PICARD_SORTSAM(BWA_MEM.out.sam)

  // Mark duplicates and gather alignment summary information
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  PICARD_COLLECTALIGNMENTSUMMARYMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)
  PICARD_COLLECTWGSMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)

  // 7) Collect .bam filenames in its own list
  bams = PICARD_MARKDUPLICATES.out.dedup_bam
                .collect()
  CREATE_BAMLIST(bams)

  // 8) Generate other required input files for STITCH

  if (params.do_mice) {

    CREATE_POSFILE_DO(chrs)
    stitch_inputs = CREATE_BAMLIST.out.bam_list
                                .combine(CREATE_POSFILE_DO.out.ref_files)
    RUN_STITCH_DO(stitch_inputs)
    STITCH_VCF_TO_TXT(RUN_STITCH_DO.out.stitch_vcf)
    geno_files = STITCH_VCF_TO_TXT.out.sample_genos
              .join(RUN_STITCH_DO.out.stitch_founder_genos)

  } else {

    CREATE_POSFILE(chrs)
    stitch_inputs = CREATE_BAMLIST.out.bam_list
                                .combine(CREATE_POSFILE.out.posfile)
    RUN_STITCH(stitch_inputs)
    STITCH_VCF_TO_TXT(RUN_STITCH.out.stitch_vcf)
    geno_files = STITCH_VCF_TO_TXT.out.sample_genos
              .join(RUN_STITCH.out.stitch_founder_genos)
  }
  
  STITCH_TO_QTL(geno_files)
  agg_stats = QUALITY_STATISTICS.out.quality_stats
              .join(PICARD_MARKDUPLICATES.out.dedup_metrics)
              .join(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt)
              .join(PICARD_COLLECTWGSMETRICS.out.txt)

  // may replace with multiqc
  AGGREGATE_STATS(agg_stats)
  
  markdown_template = Channel.of("${projectDir}/bin/stitch/aggregate_stats_summary.Rmd")
  align_stats = AGGREGATE_STATS.out.txt
				.collect()
  STATS_MARKDOWN(align_stats)
 
 }
