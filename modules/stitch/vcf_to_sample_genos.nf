process STITCH_VCF_TO_TXT {
  tag "$chr"

  
  errorStrategy 'ignore'
  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'
  
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "stitch.*.txt", mode:'copy'

  input:
  tuple val(chr), file(stitch_vcf)


  output:
  tuple val(chr), file("stitch.*.txt"), emit: sample_genos


  script:
  log.info "----- Converting STITCH Output VCF to Text File for Chromosome: ${chr} -----"

  """
  tabix -p vcf ${stitch_vcf}
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${stitch_vcf} > stitch.${chr}.txt
  """
}
