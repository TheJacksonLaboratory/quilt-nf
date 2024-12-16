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
  tuple file(bamlist), val(downsample_to_cov), val(chr), val(shuffle_bin_radius), val(ref_hap_dir), val(cross_type)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("quilt.*.vcf.gz"), file("quilt.*.vcf.gz.tbi"), val(ref_hap_dir), val(cross_type), emit: quilt_vcf

  script:
  log.info "----- Running QUILT on Chromosome ${chr}, ${downsample_to_cov}X, ${shuffle_bin_radius} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/run_quilt.R ${bamlist} \
      ${chr} \
      ${params.covar_file} \
      ${cross_type} \
      ${ref_hap_dir}/chr${chr}.hap.gz \
      ${ref_hap_dir}/chr${chr}.samples \
      ${ref_hap_dir}/chr${chr}.legend.gz \
      ${shuffle_bin_radius}
  """
}
