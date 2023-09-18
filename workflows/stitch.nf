#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for preparing lcWGS data for haplotype reconstruction
// using STITCH

// Import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/stitch.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {FASTQC} from "${projectDir}/modules/utility_modules/fastqc"
include {MULTIQC} from "${projectDir}/modules/utility_modules/multiqc"
include {FASTP} from "${projectDir}/modules/utility_modules/fastp"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"
include {SAMPLE_COVERAGE} from "${projectDir}/modules/samtools/calc_pileups"
include {DOWNSAMPLE_BAM} from "${projectDir}/modules/samtools/downsample_bam"
include {CREATE_BAMLIST} from "${projectDir}/modules/utility_modules/create_bamlist"
include {DO_FILTER_SANGER_SNPS} from "${projectDir}/modules/bcftools/DO_filter_sangerSNPs"
include {RUN_QUILT} from "${projectDir}/modules/quilt/run_quilt"
include {QUILT_LD_PRUNING} from "${projectDir}/modules/bcftools/quilt_LD_prune.nf"
include {QUILT_TO_QTL2} from "${projectDir}/modules/quilt/quilt_to_qtl2"
include {GENOPROBS} from "${projectDir}/modules/quilt/genoprobs"


//keep this in in case I need to revive some of the processe
//include {EXPAND_BED} from "${projectDir}/modules/utility_modules/expand_bed.nf"
//include {PILEUPS_TO_BAM} from "${projectDir}/modules/bedtools/filter_bams_to_coverage"
//include {INDEX_FILTERED_BAM} from "${projectDir}/modules/samtools/index_covered_bam"
//include {CREATE_POSFILE} from "${projectDir}/modules/bcftools/create_posfile"
//include {CREATE_POSFILE_DO} from "${projectDir}/modules/bcftools/create_posfile_DO"
//include {RUN_STITCH} from "${projectDir}/modules/stitch/run_stitch"
//include {RUN_STITCH_DO} from "${projectDir}/modules/stitch/run_stitch_DO"
//include {STITCH_VCF_TO_TXT} from "${projectDir}/modules/stitch/vcf_to_sample_genos"
//include {STITCH_TO_QTL} from "${projectDir}/modules/stitch/stitch_to_qtl2files"
//include {GENO_PROBS} from "${projectDir}/modules/stitch/genoprobs"
//include {TRIMMOMATIC_PE} from "${projectDir}/modules/utility_modules/trimmomatic"
//include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
//include {BOWTIE2} from "${projectDir}/modules/bowtie2/bowtie2"
//include {AGGREGATE_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_wgs"
//include {STATS_MARKDOWN} from "${projectDir}/modules/utility_modules/render_stats_markdown"
//include {GATK_HAPLOTYPECALLER_INTERVAL} from "${projectDir}/modules/gatk/gatk_haplotypecaller_interval.nf"
//include {COMBINE_GVCF} from "${projectDir}/modules/gatk/combine_gvcfs.nf"
//include {GENOTYPE_COMBINED_GVCF} from "${projectDir}/modules/gatk/genotype_combined_gvcfs.nf"
//include {GATK_VCF_TO_TXT} from "${projectDir}/modules/gatk/gatk_to_sample_genos"
//include {GATK_TO_QTL} from "${projectDir}/modules/gatk/gatk_to_qtl2"
//include {WRITE_QTL2_FILES} from "${projectDir}/modules/gatk/write_qtl2files"
//include {GENO_PROBS} from "${projectDir}/modules/gatk/genoprobs"
//include {MAKE_B6_VARIANTS} from "${projectDir}/modules/quilt/make_B6_sanger_variants"
//include {MAKE_QUILT_REFERENCE_FILES} from "${projectDir}/modules/quilt/make_haplegendsample"
//include {MAKE_QUILT_MAP} from "${projectDir}/modules/quilt/make_quilt_map"





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

chrs = Channel.of(1..19,"X")

// main workflow
workflow QUILT {

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
  //TRIMMOMATIC_PE(read_ch)

  // Calculate quality statistics for sequencing
  FASTP(read_ch)
  //QUALITY_STATISTICS(read_ch)
  
  // Run fastqc on adapter trimmed reads
  FASTQC(FASTP.out.fastp_filtered)

  // Run multiqc
  fastqc_reports = FASTQC.out.to_multiqc.flatten().collect()
  MULTIQC(fastqc_reports)

  // Generate read groups
  READ_GROUPS(FASTP.out.fastp_filtered, "gatk")
  mapping = FASTP.out.fastp_filtered.join(READ_GROUPS.out.read_groups)

  // Alignment
  BWA_MEM(mapping)
  // BOWTIE2(mapping)

  // Sort SAM files
  PICARD_SORTSAM(BWA_MEM.out.sam)

  // Mark duplicates 
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  data = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)

  // Calculate pileups
  PICARD_COLLECTALIGNMENTSUMMARYMETRICS(data)
  PICARD_COLLECTWGSMETRICS(data)
  SAMPLE_COVERAGE(data)
  
  // Downsample bams to specified coverage if the full coverage allows
  coverageFilesChannel = SAMPLE_COVERAGE.out.depth_out.map { 
	tuple -> [tuple[0], tuple[1].splitText()[0].replaceAll("\\n", "").toFloat()] 
  }

  // downsample bam files
  DOWNSAMPLE_BAM(coverageFilesChannel.join(SAMPLE_COVERAGE.out.bam_out))

  // Collect downsampled .bam filenames in its own list
  bams = DOWNSAMPLE_BAM.out.downsampled_bam.collect()
  CREATE_BAMLIST(bams)

  // (No downsampling of bams)
  // Collect .bam filenames in its own list
  //bams = PICARD_MARKDUPLICATES.out.dedup_bam.collect()
  //CREATE_BAMLIST(bams)
  
  // Run QUILT
  quilt_inputs = CREATE_BAMLIST.out.bam_list.combine(chrs)
  RUN_QUILT(quilt_inputs)

  // Perform LD pruning on QUILT output
  //QUILT_LD_PRUNING(RUN_QUILT.out.quilt_vcf)

  // Convert QUILT outputs to qtl2 files
  quilt_for_qtl2 = RUN_QUILT.out.quilt_vcf
  QUILT_TO_QTL2(quilt_for_qtl2)

  // Reconstruct haplotypes with qtl2
  GENOPROBS(QUILT_TO_QTL2.out.qtl2files)
  
 }

