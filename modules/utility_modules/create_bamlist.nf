process CREATE_BAMLIST {

  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'quay.io-jaxcompsci-rstudio-4.2.0'

  publishDir "${params.sample_folder}/bams", pattern: "STITCH_bamlist.txt", mode:'copy'

  input:
  val(bams)

  output:
  path('STITCH_bamlist.txt'), emit: bam_list

  script:
  log.info "----- Create List of .bam Files for STITCH -----"

  """
  echo ${bams} > STITCH_bamlist.txt
  Rscript --vanilla ${projectDir}/bin/shared/create_bamlist.R
  """
}
