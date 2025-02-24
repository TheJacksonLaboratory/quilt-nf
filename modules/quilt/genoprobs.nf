process GENOPROBS {
  tag "$chr, $downsample_to_cov"

  cpus 2
  memory {200.GB * task.attempt}
  time {12.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  input:
  tuple val(chr), val(downsample_to_cov), val(start), val(stop), val(shuffle_bin_radius), file(founder_geno), file(sample_genos), file(pmap), file(gmap), file(covar), file(pheno)

  output:
  tuple val(chr), val(downsample_to_cov), val(start), val(stop), val(shuffle_bin_radius), file("*36_state_probs.RData"), file("*_cross.RData"), emit: geno_probs_out, optional: true

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/genoprobs.R ${chr} \
	${sample_genos} \
	${founder_geno} \
	${pmap} \
	${gmap} \
	${covar} \
	${params.cross_type}
  """
}
