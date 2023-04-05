process GENO_PROBS {

  ag "$sampleID"

  cpus 8
  memory 400.GB
  time '24:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  input:
  tuple val(chr), val(sampleID), file(sample_geno), file(allele_codes), file(pmap), file(gmap), file(foundergeno)

  output:
  tuple val(chr), file("*_crossinfo.csv") ,file("*36_state_probs.RData"), file("*8_state_probs.RData"), file("*_control_file.json"), emit: geno_probs_out

  script:
  log.info "----- Reconstructing Sample Haplotypes for ${sampleID} -----"

  """
  Rscript --vanilla ${projectDir}/bin/gatk/genoprobs.R ${sampleID} ${params.covar_file} ${params.sample_folder}
  mv *36_state_probs.RData ${sampleID}_36_state_probs.RData
  mv *8_state_probs.RData ${sampleID}_8_state_probs.RData
  """
}
