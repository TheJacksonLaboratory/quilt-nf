process BOWTIE2 {
  tag "$sampleID"

  cpus 10
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'docker://biocontainers/bowtie2:v2.4.1_cv1'

  //publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.sam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(fq_reads), file(read_groups)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam

  script:
  log.info "----- Bowtie2 Alignment Running on: ${sampleID} -----"

  //if (params.read_type == "SE"){
  //  inputfq="${fq_reads[0]}"
  //  }
  //if (params.read_type == "PE"){
  //  inputfq="${fq_reads[0]} ${fq_reads[1]}"
  // }

  """
  bowtie2 --no-unal --no-discordant --fr --end-to-end -x ${params.ref_fa_indices_bowtie} -1 ${fq_reads[0]} -2 ${fq_reads[1]} -S ${sampleID}.sam
  """
}
