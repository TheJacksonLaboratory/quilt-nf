process GATK_COMBINE_CHR_GVCF {
// Note about this module
  tag "$chrom"

  cpus = 1
  memory = 15.GB
  time = '05:30:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.sample_folder}/gvcfs", pattern: "*.g.vcf", mode:'copy'

  input:
  tuple val(chrom), tuple(sampleID), tuple(gvcf)

  output:
  tuple val(chrom), file("*.g.vcf"), emit: vcf
  //tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- Combining Chromosome ${chrom} GVCFs -----"

  sample_size = size(gvcf)

  """
  echo ${chrom}
  echo ${sample_size}
  """
}
