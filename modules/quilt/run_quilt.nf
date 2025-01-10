process RUN_QUILT {
  tag "$chr, $downsample_to_cov"
  
  memory 300.GB
  time {16.hour * task.attempt}
  cpus 1
  errorStrategy 'retry' 
  maxRetries 1

  
  container 'sjwidmay-stitch_nf-QUILT'

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
      ${params.covar_file} \
      ${params.cross_type} \
      ${projectDir}/reference_data/${params.cross_type}/chr${chr}.hap.gz \
      ${projectDir}/reference_data/${params.cross_type}/chr${chr}.samples \
      ${projectDir}/reference_data/${params.cross_type}/chr${chr}.legend.gz \
      ${shuffle_bin_radius}
  """
}
