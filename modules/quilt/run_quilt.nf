process RUN_QUILT {
  tag "$chr, $downsample_to_cov"
  
  memory 200.GB
  time {16.hour * task.attempt}
  cpus 1
  errorStrategy 'retry' 
  maxRetries 2

  
  container 'docker://sjwidmay/stitch_nf:QUILT'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/quilt_vcfs", pattern:"*", mode:'copy'
  
  input:
  tuple file(bamlist), val(downsample_to_cov), val(chr), val(shuffle_bin_radius)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("quilt.*.vcf.gz"), file("quilt.*.vcf.gz.tbi"), emit: quilt_vcf

  script:
  log.info "----- Running QUILT on Chromosome ${chr}, ${downsample_to_cov}X, ${shuffle_bin_radius} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/run_quilt.R ${bamlist} \
      ${chr} \
      ${params.ref_file_dir}/chr${chr}.hap.gz \
      ${params.ref_file_dir}/chr${chr}.samples \
      ${params.ref_file_dir}/chr${chr}.legend.gz \
      ${params.covar_file} \
      ${shuffle_bin_radius}
  """
}
