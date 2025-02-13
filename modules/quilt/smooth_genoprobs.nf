process SMOOTH_GENOPROBS {

  cpus 1
  memory {200.GB * task.attempt}
  time {1.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/geno_probs", pattern:"*_smooth.rds", mode:'copy', overwrite: true
  
  input:
  tuple val(chrs), val(downsample_to_cov), val(shuffle_bin_radius), file(cross), file(geno_probs), file(allele_probs)

  output:
  tuple val(chrs), val(downsample_to_cov), val(shuffle_bin_radius), file("complete_cross_smooth.rds"), file("complete_genoprobs_smooth.rds"), file("complete_alleleprobs_smooth.rds"), emit: smooth_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/smooth_genoprobs.R ${geno_probs} \
  ${cross} \
  ${downsample_to_cov} \
  ${shuffle_bin_radius} \
  ${params.hap_block_cor_thresh} \
  ${projectDir}/bin/quilt/get_hap_blocks.R
  """
}