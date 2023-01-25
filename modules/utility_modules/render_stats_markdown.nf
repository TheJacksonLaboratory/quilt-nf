process WGS_STATS_MARKDOWN {

  cpus = 1
  time = '00:30:00'

  container 'docker://sjwidmay/stitch_nf:stats_markdown'

  publishDir "${params.sample_folder}/alignment_stats", pattern:"*.html", mode:'copy'
  publishDir "${params.sample_folder}/alignment_stats", pattern:"*.Rmd", mode:'copy'

  input:
  val(txts)

  output:
  tuple file('*.html'), file('*.Rmd') emit: markdown
  
  script:
  log.info "----- Rendering Alignment Statistics Summary -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/render_markdown.R ${projectDir}/bin/stitch/aggregate_stats_summary.Rmd
  """
}