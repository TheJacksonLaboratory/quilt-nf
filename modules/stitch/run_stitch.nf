process RUN_STITCH {
  tag "$chr"

  label "STITCH"
  
  cpus 1
  memory 50.GB
  time '00:30:00'
  
  input:
  tuple val(chr), file(posfile)

  output:
  tuple val(chr), file("*.vcf"), emit: stitch_output 

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch.R ${params.sample_folder} ${posfile} ${params.nFounders} ${chr}
  """
}