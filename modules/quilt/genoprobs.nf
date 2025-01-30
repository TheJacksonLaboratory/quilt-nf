process GENOPROBS {
  tag "$chr, $downsample_to_cov"

  cpus 1
  memory {400.GB * task.attempt}
  time {11.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1


  container 'docker://sjwidmay/lcgbs_hr:latest'

  // publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/geno_probs", pattern:"*.RData", mode:'copy'
  
  input:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file(founder_geno), file(sample_genos), file(pmap), file(gmap), file(covar), file(pheno)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("*36_state_probs.RData"), file("*_cross.RData"), emit: geno_probs_out

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
