process MAKE_B6_VARIANTS {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'
  stageOutMode 'move'

  container 'docker://sjwidmay/lcgbs_hr:variantannotation'

  input:
  tuple val(chr), file(vcf), file(vcf_index)

  output:
  tuple val(chr), file("*_do_snps.vcf.bgz"), file("*_do_snps.vcf.bgz.tbi"), file("C57BL_6J.tab"), emit: filtered_sanger_vcfs

  script:
  log.info "----- Make B6 Sanger SNPs for Chromosome ${chr} -----"

  """
  Rscript ${projectDir}/bin/stitch/filter_sanger_file.R ${vcf}
  """
}
