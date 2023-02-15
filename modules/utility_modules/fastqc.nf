process FASTQC {
    
  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '01:00:00'

  container 'quay.io-biocontainers-fastqc-0.11.9--hdfd78af_1'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple file("*_fastqc.html"), file("*_fastqc.zip"), file("*_fastqc.gz"), emit: to_multiqc

  script:
  log.info "----- FASTQC Running on Sample: ${sampleID} -----"

  """
  fastqc ${fq_reads[0]} ${fq_reads[1]} 
  """
}