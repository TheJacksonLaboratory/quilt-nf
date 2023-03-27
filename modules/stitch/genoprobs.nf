process GENO_PROBS {
  tag "$chr"

  cpus 8
  memory 400.GB
  time '24:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/geno_probs", pattern: "*.RData", mode:'copy'
  publishDir "${params.sample_folder}/qtl2files", pattern: "*_crossinfo.csv", mode:'copy'
  publishDir "${params.sample_folder}/qtl2files", pattern: "*.json", mode:'copy'

  input:
  tuple val(chr), file(sample_genos), file(allele_codes), file(pmap), file(gmap), file(founder_geno)

  output:
  tuple val(chr), file("*_crossinfo.csv") ,file("*36_state_probs.RData"), file("*8_state_probs.RData"), file("*_control_file.json"), emit: geno_probs_out

  script:
  log.info "----- Reconstructing Sample Haplotypes for Chromosome: ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/genoprobs.R ${chr} ${sample_genos} ${allele_codes} ${pmap} ${gmap} ${founder_geno} ${params.covar_file}
  mv *36_state_probs.RData chr${chr}_36_state_probs.RData
  mv *8_state_probs.RData chr${chr}_8_state_probs.RData
  """
}
