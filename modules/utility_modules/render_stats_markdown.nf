process STATS_MARKDOWN {
  cpus = 8
  memory 400.GB
  time = '02:30:00'

  container 'docker://sjwidmay/stitch_nf:stats_markdown'

  publishDir "${params.sample_folder}/alignment_stats", pattern:"*.html", mode:'copy'
  publishDir "${params.sample_folder}/alignment_stats", pattern:"*.Rmd", mode:'copy'

  input:
  tuple val(txts), file(Rmd)

  output:
  tuple file('*.html'), file('*.Rmd'), emit: markdown
  
  script:
  log.info "----- Rendering Alignment Statistics Summary -----"

  """
  cat ${Rmd} | sed 's+ALIGN_SUMMARIES+"c(${txts})"+g' > aggregate_stats_summary_working.Rmd
  Rscript --vanilla ${projectDir}/bin/stitch/render_markdown.R aggregate_stats_summary_working.Rmd
  """

}
