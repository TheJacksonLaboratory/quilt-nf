process MERGE_CHROMS {

  cpus 1
  memory 10.GB
  time 1.hour

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  input: 
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file(genoprobs), file(maps), file(original_crosses)

  output:
  tuple val(downsample_to_cov), val(shuffle_bin_radius), file("*_genoprobs.rds"), file("*_pmap.rds"), emit: chr_merged_probs

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/merge_chrom_probs.R ${chr}
  """

}