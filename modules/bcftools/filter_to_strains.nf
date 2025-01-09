process FILTER_TO_STRAINS {

  memory 50.GB
  time '1:00:00'
  
  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  // publishDir "${projectDir}/reference_data/", pattern:"GRCm39.fa", mode:'copy', overwrite: false
  // publishDir "${projectDir}/reference_data/", pattern:"*.vcf.gz", mode:'copy', overwrite: false
  
  input:
  tuple file(ref_genome), file(sanger_snps), val(strains), val(chr)

  output:
  tuple file(ref_genome), file(sanger_snps), file("mgp_REL2021_snps.vcf.gz.tbi"), val(strains), val(chr), file("*hom_seg_snps_indels.vcf.gz"), emit: filtered_sanger_snps
  
  script:
  
  """
  # Filter reference SNPs to SNPs on the chromosome
  bcftools index -t ${sanger_snps}
  bcftools filter --regions ${chr} \
      --include 'INFO/INDEL=0 && FILTER="PASS" && TYPE="snp"' \
      --output-type z \
      --output snps_indels_chr${chr}.vcf.gz \
      ${sanger_snps}
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