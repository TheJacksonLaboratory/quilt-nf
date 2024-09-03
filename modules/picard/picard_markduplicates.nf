process PICARD_MARKDUPLICATES {
  tag "$sampleID"

  cpus 1
  memory 200.GB
  time '12:00:00'

  container 'quay.io-biocontainers-picard-2.26.10--hdfd78af_0'

  // save if mouse and wes or save if keep intermediate
  //publishDir "${params.sample_folder}/bams"
  //publishDir "${params.sample_folder}/bams", pattern: "*.bam", mode:'copy'
  //publishDir "${params.sample_folder}/bams", pattern: "*.bai", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_dedup.bam"), emit: dedup_bam
  tuple val(sampleID), file("*_dedup.bai"), emit: dedup_bai
  tuple val(sampleID), file("*.txt"), emit: dedup_metrics

  script:
  log.info "----- Picard SortSam Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  picard -Xmx${my_mem}G MarkDuplicates \\
  I=${bam[0]} \\
  O=${sampleID}.sorted.marked4_dedup.bam \\
  M=${sampleID}.sorted.metrics.txt \\
  REMOVE_DUPLICATES=true \\
  CREATE_INDEX=true \\
  VALIDATION_STRINGENCY=LENIENT \\
  > ${sampleID}.picard.log 2>&1  
  """
}
