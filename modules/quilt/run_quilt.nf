process RUN_QUILT {
  tag "$chr"
  
  // Chr Y doesn't work for some reason
  // errorStrategy 'ignore'
  memory 501.GB
  time '05:00:00'
  
  container 'docker://sjwidmay/stitch_nf:QUILT'

publishDir "${params.pubdir}/${params.run_name}/quilt_vcfs", pattern:"*", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr), file(hapfile), file(legendfile), file(samples)

  output:
  tuple val(chr), file("quilt.*.vcf.gz"), file("quilt.*.vcf.gz.tbi"), emit: quilt_vcf

  script:
  log.info "----- Running QUILT on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/run_quilt_DO.R ${bamlist} ${chr} ${hapfile} ${samples} ${legendfile}
  """
}
