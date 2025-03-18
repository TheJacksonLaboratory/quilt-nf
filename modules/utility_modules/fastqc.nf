process FASTQC {
    
  tag "$sampleID"

  cpus 2
  memory 50.GB
  time '01:00:00'

  container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'

  input:
  tuple val(sampleID), file(fq_reads), file(fastp_reports)

  output:
  tuple file("*_fastqc.html"), file("*_fastqc.zip"), emit: to_multiqc

  script:

  """
  fastqc ${fq_reads}
  """
}
