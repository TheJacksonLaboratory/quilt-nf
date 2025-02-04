process MPILEUP {
  tag "$sampleID"

  cpus 1
  memory 75.GB
  time '03:00:00'

  container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

  publishDir "${params.pubdir}/${params.run_name}/mpileup", pattern: "*_filtered_calls.txt", mode:'copy'
  
  input:
  tuple val(sampleID), file(bam), file(bam_bai)

  output:
  tuple val(sampleID), file("*_filtered_calls.txt"), emit: mpileup

  script:

  """
  # get pileups
  # skipping indel calling with -I flag
  bcftools mpileup -Ou -f ${params.ref_fa} ${bam} -I -o ${sampleID}_mpileup.bcf

  # call variants
  bcftools call -mv ${sampleID}_mpileup.bcf -Oz -o ${sampleID}_calls.vcf.gz
  bcftools index ${sampleID}_calls.vcf.gz

  # write the called variants
  bcftools query --print-header -e 'QUAL<=30' -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${sampleID}_calls.vcf.gz | sed 's/[[# 0-9]*]//g' | sed 's/:GT//g' > ${sampleID}_filtered_calls.txt

  """
}
