process BWA_MEM {
  tag "$sampleID"

  cpus 8
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

  input:
  tuple val(sampleID), file(fq_reads), file(report), file(read_groups)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam

  script:
  log.info "----- BWA-MEM Alignment Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    inputfq="${fq_reads[0]}"
    }
  if (params.read_type == "PE"){
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
    }

  """
  rg=\$(cat $read_groups)
  bwa mem -R \${rg} \
  -t $task.cpus ${params.mismatch_penalty} ${params.ref_fa_indices} $inputfq > ${sampleID}.sam
  """
}
