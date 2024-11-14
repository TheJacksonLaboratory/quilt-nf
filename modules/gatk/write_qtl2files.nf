process WRITE_QTL2_FILES {

  cpus 1
  memory 3.GB
  time '00:10:00'

  publishDir "${params.sample_folder}/qtl2files/${sampleID}", pattern: "*.csv", mode:'copy'

  input:
  tuple val(chr), val(sampleID), file(sample_geno), file(allele_codes), file(pmap), file(gmap), file(foundergeno)

  output:
  tuple val(chr), val(sampleID), file(sample_geno), file(allele_codes), file(pmap), file(gmap), file(foundergeno), emit: writeout

  script:
  log.info "----- Writing R/qtl2 Input Files -----"

  """
  echo ${sampleID}
  """
}
