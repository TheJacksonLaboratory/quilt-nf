process MERGE_B6_VCF {

  memory 50.GB
  time '1:00:00'
  
  container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
  
  input:
  tuple val(strains), val(final_strain_order), val(chr), file(filtered_sanger_snps), file(b6_tab), file(filtered_sanger_bgz), file(filtered_sanger_bgz_tbi)

  output:
  tuple val(chr), file("*_merged_reordered.vcf.gz"), emit: complete_vcf
  
  script:
  
  """
  # Convert C57BL/6J tab-delimited file to VCF.
  bcftools convert \
        --tsv2vcf ${b6_tab} \
        --fasta-ref ${params.ref_fa} \
        --samples C57BL_6J \
        --output-type z \
        --output chr${chr}_C57BL_6J.vcf.gz

  # Index the C57BL/6J file.
  bcftools index chr${chr}_C57BL_6J.vcf.gz

  # Merge C57BL/6J into filtered Sanger VCF.
  bcftools merge \
        --output-type z \
        --output chr${chr}_merged.vcf.gz \
        chr${chr}_C57BL_6J.vcf.gz \
        ${filtered_sanger_bgz}

  # Reorder merged VCF.
  bcftools view \
        --samples ${final_strain_order} \
        --output-type z \
        --output chr${chr}_merged_reordered.vcf.gz \
        chr${chr}_merged.vcf.gz
  """

  stub:

  """
  touch test_merged_reordered.vcf.gz
  """
}