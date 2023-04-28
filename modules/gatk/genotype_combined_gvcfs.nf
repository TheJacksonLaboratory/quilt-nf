process GENOTYPE_COMBINED_GVCF {
// Note about this module
  tag "$chrom"

  cpus = 4
  memory 100.GB
  time = '05:30:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.sample_folder}/gvcfs", pattern: "*_genotyped.vcf.gz", mode:'copy'

  input:
  tuple val(chrom), val(sampleID), file(combined_gvcf), file(combined_gvcf_index)

  output:
  tuple val(chrom), val(sampleID), file("*_genotyped.vcf.gz"), file("*_genotyped.vcf.gz.tbi"), emit: vcf

  script:

  log.info "----- GATK Haplotype Caller Running on Chromosome ${chrom} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" GenotypeGVCFs  \
  -R ${params.ref_fa} \
  -V ${combined_gvcf} \
  -O chr${chrom}_genotyped.vcf.gz \
  """
}
