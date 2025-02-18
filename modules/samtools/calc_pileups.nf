process SAMPLE_COVERAGE {

  cpus 1
  memory 300.GB
  time '5:00:00'
  errorStrategy 'retry' 
  maxRetries 3

  container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'

  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_GW_coverage.txt", mode:'copy'
  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_captured_site_HQ_coverage.txt", mode:'copy'
  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_summary_coverage.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai)

  output:
  tuple val(sampleID), file(bam), emit: bam_out
  tuple val(sampleID), file("*_GW_coverage.txt"), file("*_captured_site_HQ_coverage.txt"), file("*_summary_coverage.txt"), emit: depth_out

  script:

  """
  samtools depth ${bam} -a | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_GW_coverage.txt
  samtools depth ${bam} -q 20 -Q 20 | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_captured_site_HQ_coverage.txt
  samtools coverage ${bam} > ${sampleID}_summary_coverage.txt
  """
}
