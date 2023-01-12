process RUN_STITCH {
  tag "$chr"
  
  label "STITCH"
  
  memory 100.GB
  time '01:00:00'
  
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*.vcf.gz", mode:'copy'
  
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
