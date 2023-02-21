process CREATE_POSFILE_DO {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  input:
  val(chr)

  output:
  tuple val(chr), file("*_pos.txt"), file("*.hap.gz"), file("*.legend.gz"), file("*.samples"), emit: ref_files

  script:
  log.info "----- Generating Position and DO Reference File for: ${chr} -----"

  """
  bcftools view ${params.ref_vcf} --regions ${chr} -m2 -M2 -v snps | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > STITCH_${chr}_pos.txt
  
  # create a B6 vcf by changing B6/NJ alt alleles to ref
  bcftools view ${params.ref_vcf} --regions ${chr} -m2 -M2 -v snps -s C57BL_6NJ | sed 's/1\\/1/0\\/0/g' | sed 's/0\\/1/0\\/0/g' | sed 's/1\\/0/0\\/0/g' | sed 's/C57BL_6NJ/B/g' | bcftools view -Oz -o B6_${chr}.vcf.gz
  tabix -p vcf B6_${chr}.vcf.gz

  # create the other main vcf
  bcftools view ${params.ref_vcf} --regions ${chr} -m2 -M2 -v snps -s A_J,C57BL_6NJ,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ | sed 's/A_J/A/g' | sed 's/C57BL_6NJ/B/g' | sed 's/129S1_SvImJ/C/g' | sed 's/NOD_ShiLtJ/D/g' | sed 's/NZO_HlLtJ/E/g' | sed 's/CAST_EiJ/F/g' | sed 's/PWK_PhJ/G/g' | sed 's/WSB_EiJ/H/g' | bcftools view -Oz -o notB6_${chr}.vcf.gz
  tabix -p vcf notB6_${chr}.vcf.gz

  # join the two, but select the altered B6 sample field
  bcftools merge B6_${chr}.vcf.gz notB6_${chr}.vcf.gz --force-samples | bcftools view -s A,B,C,D,E,F,G,H -Oz -o merged_${chr}.vcf.gz

  # convert to hap files
  bcftools convert merged_${chr}.vcf.gz --haplegendsample merged_${chr}
  
  """
}