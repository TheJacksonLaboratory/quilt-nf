process QUILT_TO_QTL2 {
  tag "$chr, $downsample_to_cov"

  time {90.min * task.attempt}
  cpus 1
  memory {60.GB * task.attempt}
  maxRetries 1
  errorStrategy 'retry'

  container 'sjwidmay-lcgbs_hr-variantannotation'

  input:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file(sample_genos), file(sample_genos_index)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("*_founder_geno.csv"), file("*_sample_geno.csv"), file("*_pmap.csv"), file("*_gmap.csv"), file("covar.csv"), file("pheno.csv"), emit: qtl2files
  tuple val(chr), file("*_resolution_summary.csv"), emit: resolution_summary

  script:

  """
  Rscript --vanilla ${projectDir}/bin/quilt/prepare_qtl2_files.R ${projectDir}/reference_data/${params.cross_type}/chr${chr}_phased_snps.vcf.gz \
	${sample_genos} \
	${params.covar_file} \
	${params.cross_type} \
	${projectDir}/reference_data/${params.cross_type}/chr${chr}_gen_map.txt ${chr}
  """
}
