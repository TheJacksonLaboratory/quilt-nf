process CREATE_POSFILE_DO {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'
  stageOutMode 'move'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  input:
  val(chr)

  output:
  tuple val(chr), file("*_pos.txt"), file("*.hap.gz"), file("*.legend.gz"), file("*.samples"), emit: ref_files

  script:
  log.info "----- Generating Position and DO Reference File for Chromosome ${chr} -----"

  """
  # create pos file for STITCH
  bcftools view ${params.DO_vcf} --regions ${chr} -m2 -M2 -v snps | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > STITCH_${chr}_pos.txt
  
  # convert to hap files
  bcftools view ${params.DO_vcf} --regions ${chr} -m2 -M2 -v snps | bcftools convert --haplegendsample merged_${chr}
  
  """
}