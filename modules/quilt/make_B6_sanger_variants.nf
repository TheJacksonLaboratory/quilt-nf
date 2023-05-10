process MAKE_B6_VARIANTS {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'
  stageOutMode 'move'

  container 'bioconductor'

  input:
  tuple val(chr), file(vcf), file(vcf_index)

  output:
  tuple val(chr), file(vcf), file(vcf_index), file("C57BL_6J.tab"), emit: sanger_vcfs

  script:
  log.info "----- Make B6 Sanger SNPs for Chromosome ${chr} -----"

  """
  Rscript${projectDir}/bin/stitch/filter_sanger_file.R ${vcf}
  """
}