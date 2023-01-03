#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/wgs.nf"
include {getLibraryId} from "${projectDir}/bin/shared/getLibraryId.nf"
include {CONCATENATE_READS_PE} from "${projectDir}/modules/utility_modules/concatenate_reads_PE"
include {CONCATENATE_READS_SE} from "${projectDir}/modules/utility_modules/concatenate_reads_SE"
include {BWA_MEM} from "${projectDir}/modules/bwa/bwa_mem"
include {BWA_MEM_HLA} from "${projectDir}/modules/bwa/bwa_mem_hla"
include {COSMIC_ANNOTATION as COSMIC_ANNOTATION_SNP;
         COSMIC_ANNOTATION as COSMIC_ANNOTATION_INDEL} from "${projectDir}/modules/cosmic/cosmic_annotation"
include {VCF_ANNOTATE as VCF_ANNOTATE_SNP;
         VCF_ANNOTATE as VCF_ANNOTATE_INDEL} from "${projectDir}/modules/vcftools/vcf_annotate"
include {SNPEFF;
         SNPEFF as SNPEFF_SNP;
         SNPEFF as SNPEFF_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_snpeff"
include {SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_SNP;
         SNPEFF_ONEPERLINE as SNPEFF_ONEPERLINE_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpeff_oneperline"
include {SNPSIFT_EXTRACTFIELDS} from "${projectDir}/modules/snpeff_snpsift/snpsift_extractfields"
include {SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_SNP;
         SNPSIFT_DBNSFP as SNPSIFT_DBNSFP_INDEL} from "${projectDir}/modules/snpeff_snpsift/snpsift_dbnsfp"
include {AGGREGATE_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_wgs"
include {READ_GROUPS} from "${projectDir}/modules/utility_modules/read_groups"
include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
include {PICARD_SORTSAM} from "${projectDir}/modules/picard/picard_sortsam"
include {PICARD_MARKDUPLICATES} from "${projectDir}/modules/picard/picard_markduplicates"
include {PICARD_COLLECTALIGNMENTSUMMARYMETRICS} from "${projectDir}/modules/picard/picard_collectalignmentsummarymetrics"
include {PICARD_COLLECTWGSMETRICS} from "${projectDir}/modules/picard/picard_collectwgsmetrics"
include {GATK_BASERECALIBRATOR} from "${projectDir}/modules/gatk/gatk_baserecalibrator"
include {GATK_APPLYBQSR} from "${projectDir}/modules/gatk/gatk_applybqsr"
include {GATK_MERGEVCF} from "${projectDir}/modules/gatk/gatk_mergevcf"
include {GATK_MERGEVCF_LIST} from "${projectDir}/modules/gatk/gatk_mergevcf_list"
include {GATK_VARIANTANNOTATOR} from "${projectDir}/modules/gatk/gatk_variantannotator"
include {GATK_HAPLOTYPECALLER_INTERVAL} from "${projectDir}/modules/gatk/gatk_haplotypecaller_interval"
include {GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_SNP;
         GATK_SELECTVARIANTS as GATK_SELECTVARIANTS_INDEL} from "${projectDir}/modules/gatk/gatk_selectvariants"
include {GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_SNP;
         GATK_VARIANTFILTRATION as GATK_VARIANTFILTRATION_INDEL} from "${projectDir}/modules/gatk/gatk_variantfiltration"
include {MAKE_VCF_LIST} from "${projectDir}/modules/utility_modules/make_vcf_list"

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
  // Step 1: Qual_Stat
  QUALITY_STATISTICS(read_ch)

  // Step 2: Get Read Group Information
  READ_GROUPS(QUALITY_STATISTICS.out.trimmed_fastq, "gatk")

  bwa_mem_mapping = QUALITY_STATISTICS.out.trimmed_fastq.join(READ_GROUPS.out.read_groups)

  // Step 3: BWA-MEM Alignment
  if (params.gen_org=='mouse'){
    BWA_MEM(bwa_mem_mapping)
    PICARD_SORTSAM(BWA_MEM.out.sam)
  }
  if (params.gen_org=='human'){ 
  	BWA_MEM_HLA(bwa_mem_mapping)
  	PICARD_SORTSAM(BWA_MEM_HLA.out.bam)
  }

  // Step 4: Variant Preprocessing - Part 1
  PICARD_MARKDUPLICATES(PICARD_SORTSAM.out.bam)

  // If Human
  if (params.gen_org=='human'){
    GATK_BASERECALIBRATOR(PICARD_MARKDUPLICATES.out.dedup_bam)
    
    apply_bqsr = PICARD_MARKDUPLICATES.out.dedup_bam.join(GATK_BASERECALIBRATOR.out.table)

    GATK_APPLYBQSR(apply_bqsr)

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(GATK_APPLYBQSR.out.bam)
    PICARD_COLLECTWGSMETRICS(GATK_APPLYBQSR.out.bam)

    // Create a chromosome channel. HaplotypeCaller does not have multithreading so it runs faster when individual chromosomes called instead of Whole Genome
    data = GATK_APPLYBQSR.out.bam.join(GATK_APPLYBQSR.out.bai)
    
    // Read a list of contigs from parameters to provide to GATK as intervals
    // for HaplotypeCaller variant regions
    chroms = Channel
     .fromPath("${params.chrom_contigs}")
     .splitText()
     .map{it -> it.trim()}
    
    // Applies scatter intervals from above to the BQSR bam file
    chrom_channel = data.combine(chroms)
    
    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_INTERVAL(chrom_channel)
    // Gather intervals from scattered HaplotypeCaller operations into one
    // common stream for output
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_INTERVAL.out.vcf.groupTuple(),chroms.toList())
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)
    PICARD_COLLECTWGSMETRICS(PICARD_MARKDUPLICATES.out.dedup_bam)

    // create a chromosome channel. HaplotypeCaller runs faster when individual chromosomes called instead of Whole Genome
    data = PICARD_MARKDUPLICATES.out.dedup_bam.join(PICARD_MARKDUPLICATES.out.dedup_bai)
    
    // Read a list of contigs from parameters to provide to GATK as intervals
    // for HaplotypeCaller variant regions
    chroms = Channel
     .fromPath("${params.chrom_contigs}")
     .splitText()
     .map{it -> it.trim()}
    
    // Applies scatter intervals from above to the BQSR bam file
    chrom_channel = data.combine(chroms)

    // Use the Channel in HaplotypeCaller
    GATK_HAPLOTYPECALLER_INTERVAL(chrom_channel)
    // Gather intervals from scattered HaplotypeCaller operations into one
    // common stream for output
    MAKE_VCF_LIST(GATK_HAPLOTYPECALLER_INTERVAL.out.vcf.groupTuple(), chroms.toList())
    // Sort VCF within MAKE_VCF_LIST
    GATK_MERGEVCF_LIST(MAKE_VCF_LIST.out.list)
  }


  // SNP
    select_var_snp = GATK_MERGEVCF_LIST.out.vcf.join(GATK_MERGEVCF_LIST.out.idx)
    GATK_SELECTVARIANTS_SNP(select_var_snp, 'SNP')
    var_filter_snp = GATK_SELECTVARIANTS_SNP.out.vcf.join(GATK_SELECTVARIANTS_SNP.out.idx)
    GATK_VARIANTFILTRATION_SNP(var_filter_snp, 'SNP')

  // INDEL
    select_var_indel = GATK_MERGEVCF_LIST.out.vcf.join(GATK_MERGEVCF_LIST.out.idx)
    GATK_SELECTVARIANTS_INDEL(select_var_indel, 'INDEL')
    var_filter_indel = GATK_SELECTVARIANTS_INDEL.out.vcf.join(GATK_SELECTVARIANTS_INDEL.out.idx)
    GATK_VARIANTFILTRATION_INDEL(var_filter_indel, 'INDEL')

  // Cat Output to vcf-annotate* and add dbSNP annotations. 
    VCF_ANNOTATE_SNP(GATK_VARIANTFILTRATION_SNP.out.vcf, 'SNP')
    VCF_ANNOTATE_INDEL(GATK_VARIANTFILTRATION_INDEL.out.vcf, 'INDEL')
  
// Final Post-Processing Steps Differ for Human and Mouse

  // If Human
  if (params.gen_org=='human'){

    // SNP
      COSMIC_ANNOTATION_SNP(VCF_ANNOTATE_SNP.out.vcf)
      SNPEFF_SNP(COSMIC_ANNOTATION_SNP.out.vcf, 'SNP', 'vcf')
      SNPSIFT_DBNSFP_SNP(SNPEFF_SNP.out.vcf, 'SNP')
      SNPEFF_ONEPERLINE_SNP(SNPSIFT_DBNSFP_SNP.out.vcf, 'SNP')
    // INDEL
      COSMIC_ANNOTATION_INDEL(VCF_ANNOTATE_INDEL.out.vcf)
      SNPEFF_INDEL(COSMIC_ANNOTATION_INDEL.out.vcf, 'INDEL', 'vcf')
      SNPSIFT_DBNSFP_INDEL(SNPEFF_INDEL.out.vcf, 'INDEL')
      SNPEFF_ONEPERLINE_INDEL(SNPSIFT_DBNSFP_INDEL.out.vcf, 'INDEL')
      
    // Merge SNP and INDEL and Aggregate Stats
      vcf_files = SNPEFF_ONEPERLINE_SNP.out.vcf.join(SNPEFF_ONEPERLINE_INDEL.out.vcf)
      GATK_MERGEVCF(vcf_files)

      SNPSIFT_EXTRACTFIELDS(GATK_MERGEVCF.out.vcf)
  }

  // If Mouse
  if (params.gen_org=='mouse'){
    // Merge SNP and INDEL

    vcf_files = VCF_ANNOTATE_SNP.out.vcf.join(VCF_ANNOTATE_INDEL.out.vcf)

    GATK_MERGEVCF(vcf_files)

    SNPEFF(GATK_MERGEVCF.out.vcf, 'BOTH', 'gatk')

    merged_vcf_files = GATK_MERGEVCF.out.vcf.join(SNPEFF.out.vcf)

    GATK_VARIANTANNOTATOR(merged_vcf_files)

    SNPSIFT_EXTRACTFIELDS(GATK_VARIANTANNOTATOR.out.vcf)

  }

  agg_stats = QUALITY_STATISTICS.out.quality_stats.join(PICARD_MARKDUPLICATES.out.dedup_metrics).join(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.txt).join(PICARD_COLLECTWGSMETRICS.out.txt)

  // may replace with multiqc
  AGGREGATE_STATS(agg_stats)

}
