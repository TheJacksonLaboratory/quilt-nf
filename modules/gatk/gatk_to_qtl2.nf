process GATK_TO_QTL {

  cpus 1
  memory 31.GB
  time '01:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  //errorStrategy 'ignore' 

  publishDir "${params.sample_folder}/qtl2files/${sampleID}", pattern: "*.csv", mode:'copy'

  input:
  tuple val(chr), val(sampleID), file(sample_genos), file(founder_genos)

  output:
  tuple val(chr), val(sampleID), file("geno*"), file("gmap*"), file("pmap*"), file("foundergeno*"), file("allele_codes*"), emit: qtl2files

  script:
  log.info "----- Converting GATK Genotypes to R/qtl2 Input Files -----"

  """
  Rscript --vanilla ${projectDir}/bin/gatk/gatk2qtl2files.R ${chr} ${sample_genos} ${founder_genos} ${params.nFounders}
  """
}
