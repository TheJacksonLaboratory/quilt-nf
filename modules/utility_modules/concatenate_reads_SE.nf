process CONCATENATE_READS_SE {

  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '03:00:00'

  input:
  tuple val(sampleID), file(R1)

  output:
  tuple val(sampleID), file("*fastq.gz"), emit: concat_fastq

  script:
  log.info "----- Concatenate Reads Running on: ${sampleID} -----"

  """
  cat $R1 > ${sampleID}_R1.fastq.gz
  """
}
