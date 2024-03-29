
process GATK_VCF_TO_TXT {
  tag "$chrom"

  
  //errorStrategy 'ignore'
  cpus 1
  memory 51.GB
  time '00:30:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'
  
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*.txt", mode:'copy'

  input:
  tuple val(chrom), val(sampleID), file(vcf), file(vcf_index)


  output:
  tuple val(chrom), val(sampleID), file("*_gatk.txt"), file("founders_chr*.txt"), emit: sample_genos


  script:
  log.info "----- Converting GATK Output VCF to Text File for Chromosome ${chrom} -----"

  """
  bcftools view --regions ${chrom} -m2 -M2 -v snps ${params.DO_vcf} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' > founders_chr${chrom}.txt
  bcftools view --regions ${chrom} -m2 -M2 -v snps ${vcf} | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' > samples_chr${chrom}_gatk.txt
  """
}
