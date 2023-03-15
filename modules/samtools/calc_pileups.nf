process MPILEUP {

  cpus 1
  memory 30.GB
  time '1:00:00'

  container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.bam"), file("*.mpileup"), emit: mpileup

  script:
  log.info "----- Create Pileups for Sample: ${sampleID} -----"

  """
  samtools mpileup -f ${params.ref_fa} ${bam} > ${sampleID}.mpileup
  """
}