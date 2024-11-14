process RUN_STITCH_DO {
  tag "$chr"
  
  // Chr Y doesn't work for some reason
  // errorStrategy 'ignore'
  memory 400.GB
  time '48:00:00'
  
  container 'docker://sjwidmay/stitch_nf:stitch'

  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "*.vcf.gz", mode:'copy'
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "RData/EM.all.*.RData", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr), file(posfile), file(hapfile), file(legendfile), file(samples)

  output:
  tuple val(chr), file("stitch.*.vcf.gz"), emit: stitch_vcf 
  tuple val(chr), file("RData/EM.all.*.RData"), emit: stitch_founder_genos

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch_DO.R ${bamlist} ${posfile} ${params.nFounders} ${chr} ${hapfile} ${samples} ${legendfile}
  """
}
