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
  tuple val(chrom), file("*_gatk.txt"), file("founders_chr*.txt"), emit: sample_genos


  script:
  log.info "----- Converting GATK Output VCF to Text File for ${sampleID} Chromosome ${chrom} -----"

  """
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' ${vcf} > ${sampleID}_${chrom}_gatk.txt
  bcftools view --regions ${chrom} -m2 -M2 -v snps ${params.DO_vcf} > founders_chr${chrom}.vcf.gz
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' founders_chr${chrom}.vcf.gz > founders_chr${chrom}.txt
  """
}
