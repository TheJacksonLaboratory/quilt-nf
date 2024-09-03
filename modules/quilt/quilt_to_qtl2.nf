process QUILT_TO_QTL2 {
  tag "$chr, $downsample_to_cov"

  time {4.hour * task.attempt}
  cpus 1
  memory {100.GB * task.attempt}
  maxRetries 3
  errorStrategy { task.exitStatus == 138..143 ? 'retry' : 'terminate' }

  container 'sjwidmay-lcgbs_hr-variantannotation'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/qtl2files", pattern:"*", mode:'copy'

  input:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file(sample_genos), file(sample_genos_index)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("*_founder_geno.csv"), file("*_sample_geno.csv"), file("*_pmap.csv"), file("*_gmap.csv"), file("covar.csv"), file("pheno.csv"), emit: qtl2files
  tuple val(chr), file("*_resolution_summary.csv"), emit: resolution_summary

  script:
  log.info "----- Converting QUILT Genotypes to R/qtl2 Input Files for Chromosome: ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/prepare_qtl2_files.R ${params.ref_file_dir}/chr${chr}_phased_snps.vcf.gz \
	${sample_genos} \
	${params.covar_file} \
	${params.cross_type} \
	${params.ref_file_dir}/chr${chr}_gen_map.txt ${chr}
  """
}
