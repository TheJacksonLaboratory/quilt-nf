process MAKE_QUILT_REFERENCE_FILES {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  input:
  tuple val(chr), file(vcf_bgz), file(vcf_bgz_index), file(tab_B6)

  output:
  tuple val(chr), file("*_do_snps.vcf.gz"), file("*_do_snps.vcf.gz.csi"), file("*.hap.gz"), file("*.legend.gz"), file("*.samples"), emit: haplegendsample

  script:
  log.info "----- Filtering DO Sanger SNPs for Chromosome ${chr} -----"

  """
  # Convert C57BL/6J tab-delimited file to VCF.
  bcftools convert \
        --tsv2vcf ${tab_B6} \
        --fasta-ref ${params.ref_fa} \
        --samples C57BL_6J \
        --output-type z \
        --output C57BL_6J.vcf.gz

  # Index the C57BL/6J file.
  bcftools index C57BL_6J.vcf.gz

  # Merge C57BL/6J into filtered Sanger VCF.
  bcftools merge \
        --output-type z \
        --output sanger_merged.vcf.gz \
        C57BL_6J.vcf.gz \
        ${vcf_bgz}


  # Make the final file phased since all of the SNPs are homozygous.
  zcat sanger_merged.vcf.gz | sed '/^##/! s/\\//\\|/g' > sanger_chr${chr}_do_snps.vcf
  bcftools view sanger_chr${chr}_do_snps.vcf | head -n100
  
  bgzip sanger_chr${chr}_do_snps.vcf

  # Index the Sanger VCF.
  bcftools index sanger_chr${chr}_do_snps.vcf.gz

  # make haplegendsample files
  bcftools convert sanger_chr${chr}_do_snps.vcf.gz --haplegendsample sanger_chr${chr}_do

  """
}
