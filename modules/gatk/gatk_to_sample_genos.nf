process GATK_VCF_TO_TXT {
  tag "$chr"

  
  //errorStrategy 'ignore'
  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'
  
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*_gatk.txt", mode:'copy'

  input:
  tuple val(chrom), val(sampleID), file(vcf)


  output:
  tuple val(chr), file("*_gatk.txt"), emit: sample_gneos


  script:
  log.info "----- Converting GATK Output VCF to Text File for ${sampleID} Chromosome ${chrom} -----"

  """
  tabix -p vcf ${vcf}
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${vcf} > ${sampleID}_${chrom}_gatk.txt
  """
}
