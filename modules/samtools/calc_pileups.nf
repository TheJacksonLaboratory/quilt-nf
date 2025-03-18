process SAMPLE_COVERAGE {

  cpus 1
  memory 10.GB
  time '3:00:00'
  errorStrategy 'retry' 
  maxRetries 3

  container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_GW_coverage.txt", mode:'copy', overwrite: true
  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_summary_coverage.txt", mode:'copy', overwrite: true
  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_chrX_coverage.txt", mode:'copy', overwrite: true

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file(bam), emit: bam_out
  tuple val(sampleID), file("*_GW_coverage.txt"), file("*_summary_coverage.txt"), emit: depth_out
  tuple file("*_GW_coverage.txt"), file("*_chrX_coverage.txt"), emit: chrX_depth_out

  script:

  """
  samtools depth ${bam} -a | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_GW_coverage.txt
  samtools depth ${bam} -a -r X | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_chrX_coverage.txt
  samtools coverage ${bam} > ${sampleID}_summary_coverage.txt
  """
}
