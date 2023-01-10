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
  log.info "----- Generating Position File for: ${chr} -----"

  """
  bcftools view ${params.ref_vcf} --regions ${chr} | bcftools query -f '%CHROM %POS %REF %ALT\n' | awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}'  > STITCH_${chr}_pos.txt

  """
}
