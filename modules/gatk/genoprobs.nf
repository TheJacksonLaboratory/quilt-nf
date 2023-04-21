process GENO_PROBS {

  cpus 8
  memory 501.GB
  time '24:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  input:
  val(chr)  
  output:
  val(chr)  

  script:
  log.info "----- Reconstructing Sample Haplotypes -----"

  """
  Rscript --vanilla ${projectDir}/bin/gatk/genoprobs.R ${params.covar_file} ${params.sample_folder}
  """
}
