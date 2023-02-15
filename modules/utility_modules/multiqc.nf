process MULTIQC {

  cpus 1
  memory 30.GB
  time '01:00:00'

  container 'quay.io-biocontainers-multiqc-1.12--pyhdfd78af_0'

  input:
  tuple file(fastqc_html), file(fastqc_zip), file(fastqc_gz)

  output:
  tuple file("*_fastqc.html"), emit: multiqc_report

  script:
  log.info "----- FASTQC Running on Sample: ${sampleID} -----"

  """
  multiqc .
  """
}