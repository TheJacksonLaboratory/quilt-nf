process CONCATENATE_GENOPROBS {

  cpus 2
  memory {400.GB * task.attempt}
  time {3.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.pubdir}/${params.run_name}/${shuffle_bin_radius}/geno_probs", pattern:"*.rds", mode:'copy'
  
  input:
  tuple val(chrs), val(downsample_to_cov), val(shuffle_bin_radius), file(genoprobs), file(crosses), file(interp_gridfile)

  output:
  tuple val(downsample_to_cov), val(shuffle_bin_radius), file("complete_pmap.rds"), file("complete_genoprobs.rds"), file("complete_alleleprobs.rds"), file("interp_250k_alleleprobs.rds"), file("grid_pmap.rds"), emit: concat_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/concatenate_genoprobs.R ${params.cross_type} ${interp_gridfile} ${projectDir}/bin/quilt/interpolate_genoprobs.R
  """
}
