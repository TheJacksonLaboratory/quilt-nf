process BWA_MEM_DDRADSEQ {
  tag "$sampleID"

  cpus 8
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1'

  //publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'bwa_mem' }", pattern: "*.sam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(fq_1), file(fq_2), file(read_groups)

  // for runs without clone filtering
  //input:
  //tuple val(sampleID), file(fq_files), file(read_groups)

  output:
  tuple val(sampleID), file("*.sam"), emit: sam

  script:
  log.info "----- BWA-MEM Alignment Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    inputfq="${fq_1}"
    }
  if (params.read_type == "PE"){
    inputfq="${fq_1} ${fq_2}"
    }

  """
  rg=\$(cat $read_groups)
  bwa mem -R \${rg} \
  -t $task.cpus ${params.mismatch_penalty} ${params.ref_fa_indices} $inputfq > ${sampleID}.sam
  """
}
