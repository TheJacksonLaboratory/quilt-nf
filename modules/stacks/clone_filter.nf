process CLONE_FILTER {

  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '02:00:00'

  container 'docker://sjwidmay/stitch_nf:stacks'

  //publishDir "${params.sample_folder}/fastp", pattern:"*_fastp_report.html", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*1.1.fq.gz*"), file("*2.2.fq.gz"), emit: clone_filtered

  script:
  log.info "----- Stacks Clone Filtering: ${sampleID} -----"


  """
  /stacks-2.64/clone_filter -1 ${fq_reads[0]} \\
                            -2 ${fq_reads[1]} \\
                            -i gzfastq \\
                            -y gzfastq \\
                            -o ${outDir} \\
                            -D \\
                            --oligo_len_1 8 \\
                            --inline_null
  """
}