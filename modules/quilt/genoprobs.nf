process GENOPROBS {
  tag "$chr, $downsample_to_cov"

  cpus 2
  memory {200.GB * task.attempt}
  time {12.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://sjwidmay/lcgbs_hr:latest'

  publishDir "${params.pubdir}/${params.run_name}/${shuffle_bin_radius}/geno_probs", pattern:"*.RData", mode:'copy'
  
  input:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), path(founder_geno), path(sample_genos), path(pmap), path(gmap), path(covar), path(pheno)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), path("*36_state_probs.RData"), path("*_cross.RData"), emit: geno_probs_out

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
