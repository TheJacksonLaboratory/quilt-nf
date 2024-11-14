process PICARD_COLLECTWGSMETRICS {
  tag "$sampleID"

  cpus = 1
  memory = 50.GB
  time = '08:00:00'

  container 'broadinstitute-gatk-4.2.4.1'

  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern: "*.txt", mode:'copy'
  publishDir "${params.pubdir}/${params.run_name}/bams", pattern: "*.bam", mode:'copy'
  publishDir "${params.pubdir}/${params.run_name}/bams", pattern: "*.bai", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bam_bai)

  output:
  tuple val(sampleID), file("*.txt"), emit: txt

  script:
  log.info "----- Collect WGS Metrics Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G" CollectWgsMetrics \
  --INPUT ${bam} \
  --OUTPUT ${sampleID}_CollectWgsMetrics.txt \
  --REFERENCE_SEQUENCE ${params.ref_fa} \
  --VALIDATION_STRINGENCY LENIENT
  """
}
