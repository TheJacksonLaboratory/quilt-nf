process MAKE_QUILT_MAP {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'

  container 'docker://sjwidmay/lcgbs_hr:variantannotation'

  input:
  tuple val(chr), file(filtered_vcf), file(filtered_vcf_ind)

  output:
  tuple val(chr), file(filtered_vcf), file(filtered_vcf_ind), file("*_gen_map.txt"), emit: quilt_map
  
  script:
  log.info "----- Making QUILT Genetic Map File for Chromosome ${chr} -----"

  """
  Rscript ${projectDir}/bin/quilt/make_quilt_map_file.R ${filtered_vcf} ${chr}
  """
}
