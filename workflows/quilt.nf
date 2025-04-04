#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for preparing lcWGS data for haplotype reconstruction using QUILT

// Import modules
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
include {SAMPLE_COVERAGE} from "${projectDir}/modules/samtools/calc_pileups"
include {SEX_CHECK} from "${projectDir}/modules/utility_modules/sex_check"
include {DOWNSAMPLE_BAM} from "${projectDir}/modules/samtools/downsample_bam"
include {CREATE_BAMLIST} from "${projectDir}/modules/utility_modules/create_bamlist"
include {RUN_QUILT} from "${projectDir}/modules/quilt/run_quilt"
include {QUILT_TO_QTL2} from "${projectDir}/modules/quilt/quilt_to_qtl2"
include {GENOPROBS} from "${projectDir}/modules/quilt/genoprobs"
include {CONCATENATE_GENOPROBS} from "${projectDir}/modules/quilt/concatenate_genoprobs"

//include {SMOOTH_GENOPROBS} from "${projectDir}/modules/quilt/smooth_genoprobs"

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
        //FASTQC(FASTP.out.fastp_filtered)
        //fastqc_reports = FASTQC.out.to_multiqc.flatten().collect()

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
  to_multiqc = align_summaries
                    .mix(wgs_summaries)
                    .flatten()
                    .collect()
  MULTIQC(to_multiqc)
      
  // Calculate pileups
  SAMPLE_COVERAGE(data)
  coverage_files_channel=SAMPLE_COVERAGE.out.chrX_depth_out.collect()
  SEX_CHECK(coverage_files_channel)

  if (params.align_only == false){

  if (params.downsample){
    // Extract samtools coverage estimate    
    coverageFilesChannel = SAMPLE_COVERAGE.out.depth_out.map { 
        tuple -> [tuple[0], tuple[1].splitText()[0].replaceAll("\\n", "").toFloat()] 
    }
    
    // Use full coverage estimate and specified downsampling levels to subset the bam file
    downsampleChannel = Channel.fromPath("${params.downsample_to_cov}").splitCsv()
    downsampling_bams = coverageFilesChannel.join(SAMPLE_COVERAGE.out.bam_out).combine(downsampleChannel)
    DOWNSAMPLE_BAM(downsampling_bams)
    bams = DOWNSAMPLE_BAM.out.downsampled_bam.groupTuple(by: 1)
  } else {
    bams = PICARD_MARKDUPLICATES.out.dedup_bam.map {
      tuple -> tuple[1]
    }.collect().map{
      bamlist -> [bamlist, "no_downsample"]
    }
  }
    
  // Collect downsampled .bam filenames in its own list
  CREATE_BAMLIST(bams)
  
  // bin shuffle radius channel import
  binShuffleChannel = Channel.fromPath("${params.bin_shuffling_file}").splitCsv()
  chrChunks = Channel.fromPath("${projectDir}/reference_data/${params.cross_type}/chromosome_chunks.csv")
                    .splitCsv(header: true)
                    .map {row -> 
                            [ chr         = row.chr,
                              chunk_start = row.start,
                              chunk_stop  = row.stop] }
                    .map {it -> [ it[0].toString(), it[1], it[2] ]}
  
  // Run QUILT
  quilt_inputs = CREATE_BAMLIST.out.bam_list.combine(chrChunks).combine(binShuffleChannel).combine(SEX_CHECK.out.sex_checked_covar)
  RUN_QUILT(quilt_inputs)
  
  // Convert QUILT outputs to qtl2 files by chromosome
  quilt_for_qtl2 = RUN_QUILT.out.quilt_vcf.groupTuple()
                                          .map {tuple -> [ tuple[0], tuple[1][0], tuple[2], tuple[3], tuple[4][0], tuple[5], tuple[6], tuple[7][0] ]}
  QUILT_TO_QTL2(quilt_for_qtl2)
    
  // Reconstruct haplotypes with qtl2
  GENOPROBS(QUILT_TO_QTL2.out.qtl2files)

  // Concatenate chromosome-level genotype probs and generate whole-genome objects at 1M and 250k resolution
  interp250k_file = Channel.fromPath("${projectDir}/data/interp_0.25M_physical_grid.csv")
  collected_probs = GENOPROBS.out.geno_probs_out.groupTuple(by: [1,2]).combine(interp250k_file)
  CONCATENATE_GENOPROBS(collected_probs)

  }
}