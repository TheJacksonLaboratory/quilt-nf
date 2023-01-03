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

// import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/wgs.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"

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

// main workflow
workflow WGS {
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
  // Calculate quality statistics for sequencing
  QUALITY_STATISTICS(read_ch)

  // Generate read groups
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "gatk")
  bwa_mem_mapping = QUALITY_STATISTICS.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)

  // BWA-mem alignment
  BWA_MEM(bwa_mem_mapping)

  // Sort SAM files
  PICARD_SORTSAM(BWA_MEM.out.sam)
  // Mark duplicates and gather alignment summary information
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)
  PICARD_COLLECTALIGNMENTSUMMARYMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)
  PICARD_COLLECTWGSMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)

  // 7) Collect .bam filenames in its own list

  
  // 8) Generate other required input files for STITCH
  STITCH_INPUT(...)
  // 9) Run STITCH 
  STITCH_RUN(...{params_nfounders},{params_downsample}, etc)
