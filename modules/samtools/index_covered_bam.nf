process INDEX_FILTERED_BAM {

  cpus 1
  memory 100.GB
  time '1:00:00'
  errorStrategy 'retry' 
  maxRetries 3

  container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

  publishDir "${params.sample_folder}/bams", pattern: "*_covered.bam", mode:'copy'
  publishDir "${params.sample_folder}/bams", pattern: "*_covered.bam.bai", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file(bam), emit: covered_bam
  tuple val(sampleID), file("*_covered.bam.bai"), emit: covered_bai

  script:
  log.info "----- Index Bam for Sample: ${sampleID} -----"

  """
  samtools index -b ${bam}
  """
}
