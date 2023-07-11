process QUILT_TO_QTL2 {
  tag "$chr"

  cpus 1
  memory 15.GB
  time '01:00:00'

  container 'docker://sjwidmay/lcgbs_hr:variantannotation'

  publishDir "${params.pubdir}/${params.run_name}/qtl2files", pattern:"*", mode:'copy'

  input:
  tuple val(chr), file(sample_genos), file(sample_genos_index)

  output:
  tuple val(chr), file("*_founder_geno.csv"), file("*_sample_geno.csv"), file("*_pmap.csv"), file("*_gmap.csv"), file("covar.csv"), file("pheno.csv"), emit: qtl2files

  script:
  log.info "----- Converting QUILT Genotypes to R/qtl2 Input Files for Chromosome: ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/quilt/prepare_qtl2_files.R \
	${params.ref_file_dir}/chr${chr}_phased_snps.vcf.gz \
	${sample_genos} \
	${params.covar_file} \
	${params.cross_type} \
	${params.ref_file_dir}/chr${chr}_gen_map.txt \
	${chr}
  """
}
