process FILTER_TO_STRAINS {

  memory 50.GB
  time '1:00:00'
  
  container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

  input:
  tuple val(strains), val(final_strain_order), val(chr)

  output:
  tuple val(strains), val(final_strain_order), val(chr), file("*hom_seg_snps_indels.vcf.gz"), emit: filtered_sanger_snps
  
  script:
  
  """
  # Filter reference SNPs to SNPs on the chromosome
  bcftools filter --regions ${chr} \
      --include 'INFO/INDEL=0 && FILTER="PASS" && TYPE="snp"' \
      --output-type z \
      --output snps_indels_chr${chr}.vcf.gz \
      ${params.ref_vcf}
  bcftools index snps_indels_chr${chr}.vcf.gz

  # Filter to include only polymorphic, biallelic SNPs.
  bcftools view \
      --samples ${strains} \
       --genotype hom \
       --min-alleles 2 \
       --max-alleles 2 \
       --min-ac 2 \
       --output-type z \
       --output chr${chr}_hom_seg_snps_indels.vcf.gz \
       snps_indels_chr${chr}.vcf.gz
  """

  stub:

  """
    touch test_hom_seg_snps_indels.vcf.gz
  """
}