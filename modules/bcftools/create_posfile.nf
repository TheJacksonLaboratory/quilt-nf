process CREATE_POSFILE {
  tag "$chr"

  cpus 1
  memory 15.GB
  time '00:30:00'

  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  input:
  val(chr)

  output:
  tuple val(chr), file("*_pos.txt"), emit: posfile

  script:
  log.info "----- Generating Position File for Chromosome ${chr} -----"

  """
  # create pos file for STITCH
  bcftools view ${params.ref_vcf} --regions ${chr} -m2 -M2 -v snps | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > STITCH_${chr}_pos.txt

  """
}
