process GENO_PROBS {

  tag "$sampleID"

  cpus 8
  memory 500.GB
  time '24:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  input:
  tuple val(chr), val(sampleID), file(sample_geno), file(allele_codes), file(pmap), file(gmap), file(foundergeno)

  output:
  tuple val(sampleID), emit: geno_probs_out

  script:
  log.info "----- Reconstructing Sample Haplotypes for ${sampleID} -----"

  """
  Rscript --vanilla ${projectDir}/bin/gatk/genoprobs.R ${sampleID} ${params.covar_file} ${params.sample_folder}
  """
}
