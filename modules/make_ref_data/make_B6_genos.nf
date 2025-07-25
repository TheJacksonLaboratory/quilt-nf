process MAKE_B6_GENOS {

  memory 100.GB
  time '1:00:00'
  
  container 'docker://sjwidmay/variantannotation:latest'
  
  input:
  tuple val(strains), val(final_strain_order), val(chr), file(filtered_sanger_snps)

  output:
  tuple val(strains), val(final_strain_order), val(chr), file(filtered_sanger_snps), file("*_C57BL_6J.tab"), file("*hom_seg_snps_indels.vcf.bgz"), file("*hom_seg_snps_indels.vcf.bgz.tbi"), emit: b6_calls
  
  script:
  
  """
  # additional filtering of VCF and making B6 reference variants
  Rscript ${projectDir}/bin/make_ref_data/make_B6_variant_table.R \
        ${filtered_sanger_snps} \
        ${chr}
  """

  stub:

  """
    touch test_C57BL_6J.tab
    touch test_hom_seg_snps_indels.vcf.bgz
    touch test_hom_seg_snps_indels.vcf.bgz.tbi
  """
}