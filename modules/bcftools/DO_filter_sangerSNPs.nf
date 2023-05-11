process DO_FILTER_SANGER_SNPS {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'
  stageOutMode 'move'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  input:
  val(chr)

  output:
  tuple val(chr), file("*_do_snps.vcf.gz"), file("*_do_snps.vcf.gz.csi"), emit: sanger_vcfs

  script:
  log.info "----- Filtering DO Sanger SNPs for Chromosome ${chr} -----"

  """
  # Filter Sanger SNPs to SNPs with a PASS filter
  bcftools filter --regions ${chr} \
       --include 'INFO/INDEL=0 && FILTER="PASS" && TYPE="snp"' \
       --output-type z \
       --output s_file${chr}.vcf.gz \
       ${params.ref_vcf}

  # Filter the VCF to include only SNPs from DO founders.
  # Filter to include only polymorphic, biallelic SNPs.
  bcftools view \
       --samples A_J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ \
       --genotype hom \
       --min-alleles 2 \
       --max-alleles 2 \
       --min-ac 2 \
       --output-type z \
       --output sanger_chr${chr}_do_snps.vcf.gz \
       s_file${chr}.vcf.gz

  bcftools index sanger_chr${chr}_do_snps.vcf.gz
  """
}
