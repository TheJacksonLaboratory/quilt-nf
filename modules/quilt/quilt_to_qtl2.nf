process QUILT_TO_QTL2 {
  tag "$chr, $downsample_to_cov"

  time 60.min
  cpus 1
  memory 10.GB
  maxRetries 1
  errorStrategy 'retry'

  container 'docker://sjwidmay/bcftools_qtl2:latest'

  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/geno_probs", pattern:"covar.csv", mode:'copy', overwrite: true
  publishDir "${params.pubdir}/${params.run_name}/${downsample_to_cov}/${shuffle_bin_radius}/geno_probs", pattern:"pheno.csv", mode:'copy', overwrite: true

  input:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file(sample_genos), file(sample_genos_index)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("*_founder_geno.csv"), file("*_sample_geno.csv"), file("*_pmap.csv"), file("*_gmap.csv"), file("covar.csv"), file("pheno.csv"), emit: qtl2files
  tuple val(chr), file("*_resolution_summary.csv"), emit: resolution_summary

  script:

  """
  # make founder geno table
  bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' ${projectDir}/reference_data/${params.cross_type}/chr${chr}_phased_snps.vcf.gz | \
        sed 's/[[# 0-9]*]//g' | \
        sed 's/:GT//g' > chr${chr}_fg.txt

  # make sample geno table
  bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/HWE\\t%INFO/INFO_SCORE[\\t%GT]\\n' ${sample_genos} | \
        sed 's/[[# 0-9]*]//g' | \
        sed 's/:GT//g' > chr${chr}_sg.txt
  

  Rscript --vanilla ${projectDir}/bin/quilt/prepare_qtl2_files.R chr${chr}_fg.txt \
	  chr${chr}_sg.txt \
	  ${params.covar_file} \
	  ${params.cross_type} \
	  ${projectDir}/reference_data/${params.cross_type}/chr${chr}_gen_map.txt \
    ${chr}
  """
}
