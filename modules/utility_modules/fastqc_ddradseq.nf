process FASTQC_DDRADSEQ {
    
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container 'docker://biocontainers/fastqc:v0.11.9_cv8'

  input:
  tuple val(sampleID), file(fq_1), file(fq_2)

  output:
  tuple file("*_fastqc.html"), file("*_fastqc.zip"), emit: to_multiqc

  script:
  log.info "----- FASTQC Running on Sample: ${sampleID} -----"

  """
  fastqc ${fq_1} ${fq_2}
  """
}
