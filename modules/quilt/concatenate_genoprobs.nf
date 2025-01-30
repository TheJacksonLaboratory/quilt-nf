process CONCATENATE_GENOPROBS {

  cpus 1
  memory {50.GB * task.attempt}
  time {1.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/geno_probs", pattern:"*.rds", mode:'copy', overwrite: true
  
  input:
  tuple val(chrs), val(downsample_to_cov), val(shuffle_bin_radius), file(genoprobs), file(crosses)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("complete_cross.rds"), file("complete_genoprobs.rds"), file("complete_alleleprobs.rds"), emit: concat_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/concatenate_genoprobs.R ${params.cross_type}
  """
}
