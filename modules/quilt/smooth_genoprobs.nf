process SMOOTH_GENOPROBS {

  tag "$chr, $downsample_to_cov"

  cpus 1
  memory {320.GB * task.attempt}
  time {11.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  input:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file(geno_probs), file(cross)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("*_36_state_probs_smooth.RData"), file("*_cross_smooth.RData"), emit: smooth_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/smooth_genoprobs.R ${geno_probs} \
  ${cross} \
  ${downsample_to_cov} \
  ${shuffle_bin_radius} \
  ${params.hap_block_cor_thresh} \
  ${chr} \
  ${projectDir}/bin/quilt/get_hap_blocks.R
  """
}
