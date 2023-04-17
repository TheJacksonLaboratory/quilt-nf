process GENO_PROBS {

  tag "$sampleID"

  cpus 8
  memory 500.GB
  time '24:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  input:
  tuple val(chr), file(sample_geno), file(gmap), file(pmap), file(foundergeno), file(allele_codes)  
  output:
  tuple file("lcgbs_36_state_probs.RData"), file("lcgbs_8_state_probs.RData"), emit: geno_probs_out

  script:
  log.info "----- Reconstructing Sample Haplotypes -----"

  """
  Rscript --vanilla ${projectDir}/bin/gatk/genoprobs.R ${params.covar_file} ${params.sample_folder}
  """
}
