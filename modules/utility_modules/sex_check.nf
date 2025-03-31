process SEX_CHECK {

  cpus 1
  memory 100.GB
  time '00:30:00'

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.pubdir}/${params.run_name}/${shuffle_bin_radius}/geno_probs", pattern:"sex_check_covar.csv", mode:'copy'

  input:
  path(coverage_files)

  output:
  path('sex_check_covar.csv'), emit: sex_checked_covar

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/coverage_based_sex_check.R ${params.covar_file}
  """
}