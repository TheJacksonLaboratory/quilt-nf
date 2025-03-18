process PHASE_FOUNDER_VCF {

  memory 50.GB
  time '1:00:00'

  container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
  
  publishDir "${projectDir}/reference_data/${params.cross_type}", pattern:"*_phased_snps.vcf.gz", mode:'copy', overwrite: true
  
  input:
  tuple val(chr), file(merged_vcf)

  output:
  tuple val(chr), file("*_phased_snps.vcf.gz"), emit: phased_vcf
  
  script:
  
  """
  zcat ${merged_vcf} | sed '/^##/! s/\\//\\|/g' > chr${chr}_phased_snps.vcf
  bgzip chr${chr}_phased_snps.vcf
  """

  stub:

  """
  touch test_phased_snps.vcf
  """
}