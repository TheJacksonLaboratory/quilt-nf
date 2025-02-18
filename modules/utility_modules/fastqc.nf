process FASTQC {
    
  tag "$sampleID"

  cpus 2
  memory 50.GB
  time '01:00:00'

  container 'docker://biocontainers/fastqc:v0.11.9_cv8'

  input:
  tuple val(sampleID), file(fq_reads), file(fastp_reports)

  output:
  tuple file("*_fastqc.html"), file("*_fastqc.zip"), emit: to_multiqc

  script:

  """
  fastqc ${fq_reads[0]} ${fq_reads[1]} 
  """
}
