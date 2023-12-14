process GENOPROBS {
  tag "$chr, $downsample_to_cov"

  cpus 1
  memory {200.GB * task.attempt}
  time {12.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 3


  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/geno_probs", pattern: "*.RData", mode:'copy'
  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/geno_probs", pattern: "*_crossovers.csv", mode:'copy'  

  input:
  tuple val(chr), val(downsample_to_cov), file(founder_geno), file(sample_genos), file(pmap), file(gmap), file(covar), file(pheno)

  output:
  tuple val(chr), val(downsample_to_cov), file("*_crossovers.csv"), file("*36_state_probs.RData"), file("*8_state_probs.RData"), file("*_cross.RData"), emit: geno_probs_out

  script:
  log.info "----- Reconstructing Sample Haplotypes for Chromosome: ${chr}, ${downsample_to_cov}X -----"

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
