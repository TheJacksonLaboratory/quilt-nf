process STITCH_TO_QTL {
  tag "$chr"

  cpus 1
  memory 15.GB
  time '01:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/qtl2files", pattern: "*.csv", mode:'copy'

  input:
  tuple val(chr), file(sample_genos), file(founder_genos)

  output:
  tuple val(chr), file("geno*.csv"), file("allele_codes*.csv"), file("pmap*.csv"), file("gmap*.csv"), file("foundergeno*.csv"), emit: qtl2files

  script:
  log.info "----- Converting STITCH Genotypes to R/qtl2 Input Files for Chromosome: ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/stitch2qtl2files.R ${chr} ${sample_genos} ${params.sample_folder}/stitch_vcfs/RData/EM.all.${chr}.RData ${params.nFounders}
  """
}
