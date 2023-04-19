process GATK_TO_QTL {

  tag "$chr"  

  cpus 8
  memory 300.GB
  time '03:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  //errorStrategy 'ignore' 

  publishDir "${params.sample_folder}/qtl2files/", pattern: "*.csv", mode:'copy'

  input:
  tuple val(chr), val(sampleID), file(sample_genos), file(founder_genos)

  output:
  val(chr), emit: qtl2files

  script:
  log.info "----- Converting GATK Genotypes to R/qtl2 Input Files for Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/gatk/gatk2qtl2files.R ${chr} ${sample_genos} ${founder_genos} ${params.nFounders}
  """
}
