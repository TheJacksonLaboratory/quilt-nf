
process MULTIQC {

  cpus 1
  memory 30.GB
  time '01:00:00'

  container 'docker://ewels/multiqc:latest'
  
  publishDir "${params.pubdir}/${params.run_name}/multiqc", pattern:"*", mode:'copy'

  input:
  file('*')

  output:
  tuple file('*_multiqc_report.html'), file('multiqc_data/multiqc_data.json'), emit: multiqc_report

  script:
  log.info "----- Running MULTIQC on All Samples -----"

  """
  multiqc .
  mv multiqc_report.html ${params.run_name}_multiqc_report.html
  """
}
