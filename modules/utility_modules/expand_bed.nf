process EXPAND_BED {

  cpus 1
  memory 100.GB
  time '01:00:00'

  container 'docker://sjwidmay/lcgbs_hr:qtl2_et_al'

  publishDir "${params.sample_folder}/coverage_regions", pattern:"*_interval.bed", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bed)

  output:
  tuple val(sampleID), file(bam), file("*_interval.bed"), emit: coverage_intervals

  script:
  log.info "----- Expanding BED to Intervals for Sample ${sampleID} -----"

  """
  Rscript --vanilla ${projectDir}/bin/shared/expand_bed_for_bams.R ${bed} ${sampleID}
  """
}
