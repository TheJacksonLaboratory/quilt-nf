process RUN_QUILT {

  memory 50.GB
  time {4.hour * task.attempt}
  cpus 1
  errorStrategy 'retry' 
  maxRetries 1

  
  container 'docker://sjwidmay/quilt-nf:latest'

  input:
  tuple path(bamlist), val(downsample_to_cov), val(chr), val(start), val(stop), val(shuffle_bin_radius), path(covar_file)

  output:
  tuple val(chr), val(downsample_to_cov), val(start), val(stop), val(shuffle_bin_radius), path("quilt.*.vcf.gz"), path("quilt.*.vcf.gz.tbi"), path(covar_file), emit: quilt_vcf


  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/run_quilt.R ${bamlist} \
      ${chr} \
      ${covar_file} \
      ${params.cross_type} \
      ${projectDir}/reference_data/${params.cross_type}/chr${chr}.hap.gz \
      ${projectDir}/reference_data/${params.cross_type}/chr${chr}.samples \
      ${projectDir}/reference_data/${params.cross_type}/chr${chr}.legend.gz \
      ${shuffle_bin_radius} \
      ${start} \
      ${stop}
  """
}
