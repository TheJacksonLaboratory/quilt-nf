process CLONE_FILTER {

  tag "$sampleID"

  cpus 1
  memory 150.GB
  time '04:00:00'

  container 'docker://sjwidmay/stacks:latest'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*1.1.fq.gz*"), file("*2.2.fq.gz"), emit: clone_filtered

  script:
  
  """
  /stacks-2.64/clone_filter -1 ${fq_reads[0]} \\
                            -2 ${fq_reads[1]} \\
                            -i gzfastq \\
                            -y gzfastq \\
                            -D \\
                            --oligo_len_2 8 \\
                            --null_index
  """
}
