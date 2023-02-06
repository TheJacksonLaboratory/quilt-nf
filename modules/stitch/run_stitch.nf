process RUN_STITCH {
  tag "$chr"
  
  memory 100.GB
  time '10:00:00'
  
  container 'docker://sjwidmay/stitch_nf:latest'

  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*.vcf.gz", mode:'copy'
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "RData/EM.all.*.RData", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr), file(posfile)

  output:
  tuple val(chr), file("stitch.*.vcf.gz"), emit: stitch_vcf 
  file("RData/EM.all.*.RData"), emit: stitch_founder_genos

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch.R ${bamlist} ${posfile} ${params.nFounders} ${chr}
  """
}
