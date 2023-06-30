process RUN_QUILT {
  tag "$chr"
  
  // Chr Y doesn't work for some reason
  errorStrategy 'retry'
  maxRetries 2
  memory 645.GB
  time '06:00:00'
  cpus 1

  
  container 'docker://sjwidmay/stitch_nf:QUILT'

  publishDir "${params.pubdir}/${params.run_name}/quilt_vcfs", pattern:"*", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr)

  output:
  tuple val(chr), file("quilt.*.vcf.gz"), file("quilt.*.vcf.gz.tbi"), emit: quilt_vcf

  script:
  log.info "----- Running QUILT on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/run_quilt.R ${bamlist} \
      ${chr} \
      ${params.ref_file_dir}/chr${chr}.hap.gz \
      ${params.ref_file_dir}/chr${chr}.samples \
      ${params.ref_file_dir}/chr${chr}.legend.gz \
      ${params.covar_file}
  """
}
