process CREATE_BAMLIST {

  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'docker://rocker/rstudio:latest'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}", pattern: "bamlist.txt", mode:'copy'

  input:
  tuple val(bams), val(downsample_to_cov)

  output:
  tuple path('bamlist.txt'), val(downsample_to_cov), emit: bam_list

  script:

  """
  echo ${bams} > bamlist.txt
  Rscript --vanilla ${projectDir}/bin/quilt/create_bamlist.R
  """
}
