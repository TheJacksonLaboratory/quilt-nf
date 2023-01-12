process RUN_STITCH {
  tag "$chr"
  
  label "STITCH"
  
  cpus 16
  memory 50.GB
  time '00:30:00'
  
  errorStrategy 'ignore'

  input:
  tuple file(bamlist), val(chr), file(posfile)

  output:
  tuple val(chr), file("*.vcf"), emit: stitch_output 

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch.R ${bamlist} ${posfile} ${params.nFounders} ${chr}
  """
}
