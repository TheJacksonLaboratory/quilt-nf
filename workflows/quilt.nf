#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for preparing lcWGS data for haplotype reconstruction using QUILT

// Import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/quilt.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {FASTQC} from "${projectDir}/modules/utility_modules/fastqc"
include {MULTIQC} from "${projectDir}/modules/utility_modules/multiqc"
include {FASTP} from "${projectDir}/modules/utility_modules/fastp"
include {CLONE_FILTER} from "${projectDir}/modules/stacks/clone_filter"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"
include {MPILEUP} from "${projectDir}/modules/bcftools/mpileup"
include {SAMPLE_COVERAGE} from "${projectDir}/modules/samtools/calc_pileups"
include {DOWNSAMPLE_BAM} from "${projectDir}/modules/samtools/downsample_bam"
include {DOWNSAMPLE_BAM_ALIGN_ONLY} from "${projectDir}/modules/samtools/downsample_bam_align_only"
include {CREATE_BAMLIST} from "${projectDir}/modules/utility_modules/create_bamlist"
include {DO_FILTER_SANGER_SNPS} from "${projectDir}/modules/bcftools/DO_filter_sangerSNPs"
include {RUN_QUILT} from "${projectDir}/modules/quilt/run_quilt"
include {QUILT_TO_QTL2} from "${projectDir}/modules/quilt/quilt_to_qtl2"
include {GENOPROBS} from "${projectDir}/modules/quilt/genoprobs"


//keep this in in case I need to revive some of the processes
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
//include {QUILT_LD_PRUNING} from "${projectDir}/modules/bcftools/quilt_LD_prune.nf"
//include {FASTQC_DDRADSEQ} from "${projectDir}/modules/utility_modules/fastqc_ddradseq"
//include {BWA_MEM_DDRADSEQ} from "${projectDir}/modules/bwa/bwa_mem_ddradseq"




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

 // Concatenate Fastq files if required. 
 if (params.concat_lanes){
   if (params.read_type == 'PE'){
       CONCATENATE_READS_PE(read_ch)
       read_ch = CONCATENATE_READS_PE.out.concat_fastq
   } else if (params.read_type == 'SE'){
       CONCATENATE_READS_SE(read_ch)
       read_ch = CONCATENATE_READS_SE.out.concat_fastq
   }
 }

// Calculate quality statistics for sequencing
if (params.library_type == 'ddRADseq'){
        
        // Run Stacks clone filter
        CLONE_FILTER(read_ch)
        
        decloned_read_ch = CLONE_FILTER.out.clone_filtered.map {
  	    tuple -> [ tuple[0], [tuple[1], tuple[2]] ]
        }
        
        // remove adapters?
        FASTP(decloned_read_ch)
        
        // Run fastqc on adapter trimmed reads
        FASTQC(FASTP.out.fastp_filtered)
        fastqc_reports = FASTQC.out.to_multiqc.flatten().collect()

        // Generate read groups
        READ_GROUPS(FASTP.out.fastp_filtered, "gatk")
        mapping = FASTP.out.fastp_filtered.join(READ_GROUPS.out.read_groups)

        // Alignment
        BWA_MEM(mapping)

        // Sort SAM files
        PICARD_SORTSAM(BWA_MEM.out.sam)
        data = PICARD_SORTSAM.out.bam.join(PICARD_SORTSAM.out.bai)

    } else {

        // Trim any adapters
        FASTP(read_ch)
        
        // Run fastqc on adapter trimmed reads
        FASTQC(FASTP.out.fastp_filtered)
        fastqc_reports = FASTQC.out.to_multiqc.flatten().collect()

        // Generate read groups
        READ_GROUPS(FASTP.out.fastp_filtered, "gatk")
        mapping = FASTP.out.fastp_filtered.join(READ_GROUPS.out.read_groups)
    
        // Alignment
        BWA_MEM(mapping)

        // Sort SAM files
        PICARD_SORTSAM(BWA_MEM.out.sam)

        // Mark duplicates
        PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
        data = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
}
  // Collect reports for multiqc
  PICARD_COLLECTALIGNMENTSUMMARYMETRICS(data)
  align_summaries = PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt
  
  PICARD_COLLECTWGSMETRICS(data)
  wgs_summaries = PICARD_COLLECTWGSMETRICS.out.txt

  // Run multiqc
  to_multiqc = fastqc_reports
                    .mix(align_summaries)
                    .mix(wgs_summaries)
                    .flatten()
                    .collect()
  MULTIQC(to_multiqc)
      
  // Calculate pileups
  SAMPLE_COVERAGE(data)
  //MPILEUP(data)
  
  // Downsample bams to specified coverage if the full coverage allows
  coverageFilesChannel = SAMPLE_COVERAGE.out.depth_out.map { 
     tuple -> [tuple[0], tuple[1].splitText()[0].replaceAll("\\n", "").toFloat()] 
  }

  //downsample bam files
  downsampleChannel = Channel.fromPath("${params.downsample_to_cov}").splitCsv()
  downsampling_bams = coverageFilesChannel.join(SAMPLE_COVERAGE.out.bam_out).combine(downsampleChannel)
  DOWNSAMPLE_BAM(downsampling_bams)

  if (params.align_only == false){
    //Collect downsampled .bam filenames in its own list
    bams = DOWNSAMPLE_BAM.out.downsampled_bam.groupTuple(by: 1)
    CREATE_BAMLIST(bams)

    // For haplotype reconstruction, need to know which founder haplotypes are sourced
    if(params.cross_type == 'do'){
      ref_hap_dir = Channel.of( [params.do_ref_haps,  params.cross_type] )
    } else if(params.cross_type == 'bxd'){
      ref_hap_dir = Channel.of( [params.bxd_ref_haps, params.cross_type] )
    } else if(params.cross_type == 'cc'){
      ref_hap_dir = Channel.of( [params.do_ref_haps, params.cross_type] )
    } else {
      ref_hap_dir = Channel.of( [params.ref_hap_dir,  params.cross_type] )
    }
  
    // bin shuffle radius channel import
    binShuffleChannel = Channel.fromPath("${params.bin_shuffling_file}").splitCsv()
  
    // Run QUILT
    quilt_inputs = CREATE_BAMLIST.out.bam_list.combine(chrs)
                                              .combine(binShuffleChannel)
                                              .combine(ref_hap_dir)
    RUN_QUILT(quilt_inputs)
    // Convert QUILT outputs to qtl2 files
    quilt_for_qtl2 = RUN_QUILT.out.quilt_vcf
    QUILT_TO_QTL2(quilt_for_qtl2)
    // Reconstruct haplotypes with qtl2
    GENOPROBS(QUILT_TO_QTL2.out.qtl2files)
  }
}