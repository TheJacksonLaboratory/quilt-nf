process CREATE_BAMLIST {

  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'quay.io-jaxcompsci-rstudio-4.2.0'

  publishDir "${params.pubdir}/${params.run_name}", pattern: "bamlist.txt", mode:'copy'

  input:
  val(bams), val(downsample_to_cov)

  output:
  path('bamlist.txt'), emit: bam_list

  script:
  log.info "----- Create List of .bam Files for STITCH -----"

  """
  echo ${bams} > bamlist.txt
  Rscript --vanilla ${projectDir}/bin/stitch/create_bamlist.R
  """
}
