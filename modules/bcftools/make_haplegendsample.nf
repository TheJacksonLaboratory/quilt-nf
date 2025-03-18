process MAKE_REF_HAPLOTYPES {
  
  memory 50.GB
  time '1:00:00'

  container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

  publishDir "${projectDir}/reference_data/${params.cross_type}", pattern:"*.hap.gz", mode:'copy', overwrite: true
  publishDir "${projectDir}/reference_data/${params.cross_type}", pattern:"*.legend.gz", mode:'copy', overwrite: true
  publishDir "${projectDir}/reference_data/${params.cross_type}", pattern:"*.samples", mode:'copy', overwrite: true
  publishDir "${projectDir}/reference_data/${params.cross_type}", pattern:"*_phased_snps.vcf.gz.tbi", mode:'copy', overwrite: true
  
  input:
  tuple val(chr), file(phased_vcf)

  output:
  tuple val(chr), file("*.hap.gz"), file("*.legend.gz"), file("*.samples"), emit: haplegendsample
  tuple val(chr), file("*_phased_snps.vcf.gz.tbi"), emit: phased_vcf_tbi
  
  script:

  """
  # index vcf
  tabix -p vcf ${phased_vcf}

  # make haplegendsample files
  bcftools convert ${phased_vcf} --haplegendsample chr${chr}
  
  """

  stub:

  """
  touch test.hap.gz
  touch test.legend.gz
  touch test.samples
  """
}
