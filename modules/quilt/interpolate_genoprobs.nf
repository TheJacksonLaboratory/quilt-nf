process INTERPOLATE_GENOPROBS {

  cpus 1
  memory 10.GB
  time 1.hour

  container 'docker://sjwidmay/lcgbs_hr:latest'
  
  input: 
  tuple val(chr), val(downsample_to_cov), val(start), val(stop), val(shuffle_bin_radius), file(genoprobs), file(cross)

  output:
  tuple val(chr), val(downsample_to_cov), val(shuffle_bin_radius), file("*_pr.rds"), file("*_map.rds"), file("*_cross.rds"), emit: interpolated_probs

  script:

  """
  echo "${chr}-${start}:${stop}"

  Rscript --vanilla ${projectDir}/bin/quilt/interpolate_quilt_genoprobs.R ${chr} \
	${start} \
	${stop} \
	${genoprobs} \
	${cross} \
    ${params.gridfile} \
    ${projectDir}/bin/quilt/interpolate_genoprobs.R
  """


}