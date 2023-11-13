process RUN_QUILT {
  tag "$chr, $downsample_to_cov"
  
  memory {150.GB * task.attempt}
  time {6.hour * task.attempt}
  errorStrategy 'retry'
  maxRetries 3
  cpus 1

  
  container 'docker://sjwidmay/stitch_nf:QUILT'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/quilt_vcfs", pattern:"*", mode:'copy'
  
  input:
  tuple file(bamlist), val(downsample_to_cov), val(chr)

  output:
  tuple val(chr), val(downsample_to_cov), file("quilt.*.vcf.gz"), file("quilt.*.vcf.gz.tbi"), emit: quilt_vcf

  script:
  log.info "----- Running QUILT on Chromosome ${chr}, ${downsample_to_cov}X -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/run_quilt.R ${bamlist} \
      ${chr} \
      ${params.ref_file_dir}/chr${chr}.hap.gz \
      ${params.ref_file_dir}/chr${chr}.samples \
      ${params.ref_file_dir}/chr${chr}.legend.gz \
      ${params.covar_file}
  """
}
