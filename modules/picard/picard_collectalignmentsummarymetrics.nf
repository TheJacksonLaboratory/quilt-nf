process PICARD_COLLECTALIGNMENTSUMMARYMETRICS {

  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '03:00:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern: "*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bam_bai)

  output:
  tuple val(sampleID), file("*.txt"), emit: txt

  script:
  log.info "----- Collect Alignment Summary Metrics Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G" CollectAlignmentSummaryMetrics \
  --INPUT ${bam} \
  --OUTPUT ${sampleID}_AlignmentMetrics.txt \
  --REFERENCE_SEQUENCE ${params.ref_fa} \
  --METRIC_ACCUMULATION_LEVEL ALL_READS \
  --VALIDATION_STRINGENCY LENIENT
  """
}
