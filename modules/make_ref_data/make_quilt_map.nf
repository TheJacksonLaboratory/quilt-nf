process MAKE_QUILT_MAP {

  memory 50.GB
  time '01:00:00'

  container 'docker://sjwidmay/lcgbs_hr:variantannotation'

  publishDir "${projectDir}/reference_data/${params.cross_type}", pattern:"*_gen_map.txt", mode:'copy', overwrite: false
  
  input:
  tuple val(chr), file(phased_vcf)

  output:
  tuple val(chr), file("*_gen_map.txt"), emit: quilt_map
  
  script:

  """
  Rscript ${projectDir}/bin/quilt/make_quilt_map_file.R ${phased_vcf} ${chr}
  """

  stub:

  """
  touch test_gen_map.txt
  """
}
