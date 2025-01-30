process READ_GROUPS {
  tag "$sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  input:
  tuple val(sampleID), file(fq_reads), file(report)
  val(picard)

  output:
  tuple val(sampleID), file("*.txt"), emit: read_groups

  script:

  if (picard=="picard"){
    p='-p'
  }
  else{
    p=''
  }
  """
  python ${projectDir}/bin/shared/read_group_from_fastq.py $p -o ${sampleID}_read_group.txt ${fq_reads[0]}
  """
  }
