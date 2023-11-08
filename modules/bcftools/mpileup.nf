process MPILEUP {
  tag "$sampleID"

  cpus 1
  memory 50.GB
  time '03:00:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  publishDir "${params.pubdir}/${params.run_name}/mpileup", pattern: "*.vcf.gz", mode:'copy'
  
  input:
  tuple val(sampleID), file(bam), file(bam_bai)

  output:
  tuple val(sampleID), file("*.vcf.gz"), emit: mpileup

  script:
  log.info "----- Calling Variants with mpileup: ${sampleID} -----"

  """
  # get pileups
  bcftools mpileup -Ou -f ${params.ref_fa} ${bam} -o ${sampleID}_mpileup.bcf

  # call variants
  bcftools call -mv ${sampleID}_mpileup.bcf -Oz -o ${sampleID}_calls.vcf.gz
  bcftools index ${sampleID}_calls.vcf.gz

  """
}