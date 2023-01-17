process RUN_STITCH {
  tag "$chr"
  
  label "STITCH"
  
  memory 100.GB
  time '10:00:00'
  
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*.vcf.gz", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr), file(posfile)

  output:
  tuple val(chr), file("*.vcf.gz"), emit: stitch_output 

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch.R ${bamlist} ${posfile} ${params.nFounders} ${chr}
  """
}
