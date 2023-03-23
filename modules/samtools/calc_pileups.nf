process MPILEUP {

  cpus 1
  memory 10.GB
  time '2:00:00'
  errorStrategy 'retry' 
  maxRetries 3

  container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

  // publishDir "${params.sample_folder}/sample_coverage", pattern:"*_coverage.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file(bam), file("*.bed"), emit: bed

  script:
  log.info "----- Create Pileups for Sample: ${sampleID} -----"

  """
  samtools mpileup -f ${params.ref_fa} ${bam} > ${sampleID}.mpileup
  awk '\$4 > 5 {print \$1"\t"\$2}' ${sampleID}.mpileup > ${sampleID}.bed
  """
}
