process RUN_STITCH {
  tag "$chr"
  
  // Chr Y doesn't work for some reason
  errorStrategy 'ignore'
  memory 100.GB
  time '10:00:00'
  
  container 'docker://sjwidmay/stitch_nf:latest'

  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*.vcf.gz", mode:'copy'
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "RData/EM.all.*.RData", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr), file(posfile)

  output:
  tuple val(chr), file("stitch.*.vcf.gz"), emit: stitch_vcf 
  tuple val(chr), file("RData/EM.all.*.RData"), emit: stitch_founder_genos

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch.R ${bamlist} ${posfile} ${params.nFounders} ${chr}
  """
}
