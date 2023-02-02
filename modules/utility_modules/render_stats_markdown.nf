process STATS_MARKDOWN {
  memory 300.GB
  time = '05:00:00'

  container 'docker://sjwidmay/stitch_nf:stats_markdown'

  publishDir "${params.sample_folder}/alignment_stats", pattern:"*.html", mode:'copy'
  publishDir "${params.sample_folder}/alignment_stats", pattern:"*.Rmd", mode:'copy'

  input:
  val(txts)

  output:
  tuple file('*.html'), file('*.Rmd'), emit: markdown
  
  script:
  log.info "----- Rendering Alignment Statistics Summary -----"

  """
  echo ${txts} > summaries.txt
  cat ${projectDir}/bin/stitch/aggregate_stats_summary.Rmd > aggregate_stats_summary_working.Rmd
  Rscript ${projectDir}/bin/stitch/render_markdown.R aggregate_stats_summary_working.Rmd
  """

}
