process PILEUPS_TO_BAM {

  cpus 1
  memory 30.GB
  time '1:00:00'

  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

  publishDir "${params.sample_folder}/bams", pattern: "*_covered.bam", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bed)

  output:
  tuple val(sampleID), file("*_covered.bam"), emit: filtered_bam

  script:
  log.info "----- Filtering Alignment from ${sampleID} to Pileup Sites -----"

  """
  bedtools intersect -abam ${bam} -b ${sampleID}_interval.bed > ${sampleID}_covered.bam
  """
}
