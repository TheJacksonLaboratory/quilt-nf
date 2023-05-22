process GENOPROBS {
  tag "$chr"

  cpus 8
  memory 700.GB
  time '48:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.pubdir}/${params.run_name}/geno_probs", pattern: "*.RData", mode:'copy'
  
  input:
  tuple val(chr), file(founder_geno), file(sample_genos), file(pmap), file(gmap), file(covar), file(pheno)

  output:
  tuple val(chr), file("*36_state_probs.RData"), file("*8_state_probs.RData"), file("*_cross.RData"), emit: geno_probs_out

  script:
  log.info "----- Reconstructing Sample Haplotypes for Chromosome: ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/genoprobs.R ${chr} ${sample_genos} ${founder_geno} ${pmap} ${gmap} ${covar}
  """
}
