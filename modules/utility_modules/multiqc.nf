
process MULTIQC {

  cpus 1
  memory 30.GB
  time '01:00:00'

  container 'docker://ewels/multiqc:latest'
  
  publishDir "${params.sample_folder}/multiqc", pattern:"*", mode:'copy'
  publishDir "${params.sample_folder}/multiqc_data", pattern:"multiqc_data/*", mode:'copy'

  input:
  file('*')

  output:
  path('multiqc_report.html'), emit: multiqc_report

  script:
  log.info "----- Running MULTIQC on All Samples -----"

  """
  multiqc .
  """
}
