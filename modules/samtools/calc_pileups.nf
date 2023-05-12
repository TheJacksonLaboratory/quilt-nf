process SAMPLE_COVERAGE {

  cpus 1
  memory 300.GB
  time '5:00:00'
  errorStrategy 'retry' 
  maxRetries 3

  container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

  // publishDir "${params.sample_folder}/sample_coverage", pattern:"*_coverage.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file(bam), file("*_depth.out"), emit: depth_out

  script:
  log.info "----- Create Pileups for Sample: ${sampleID} -----"

  """
  samtools depth ${bam} -o ${sampleID}_depth.out

  # samtools mpileup -f ${params.ref_fa} ${bam} > ${sampleID}.mpileup
  # awk '\$4 > 0 {print \$1"\t"\$2}' ${sampleID}.mpileup > ${sampleID}.bed
  """
}
