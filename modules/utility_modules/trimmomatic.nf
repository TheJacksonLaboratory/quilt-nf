process TRIMMOMATIC_PE {
    
  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '24:00:00'

  container 'quay.io-biocontainers-trimmomatic-0.35--6'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'quality_stats' }", pattern: "*fastq.gz_stat", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*_paired"), emit: trimmomatic
  tuple val(sampleID), file("*_paired"), emit: to_fastqc

  script:
  log.info "----- Trimmomatic Running on: ${sampleID} -----"

  if (params.seq_method == 'lcGBS')
  """
  trimmomatic ${params.read_type} \
              ${fq_reads[0]} \
              ${fq_reads[1]} \
              ${fq_reads[0]}_paired \
              ${fq_reads[0]}_unpaired \
              ${fq_reads[1]}_paired \
              ${fq_reads[1]}_unpaired \
              ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:2:True \
              LEADING:3 TRAILING:3 MINLEN:36
  """
  else
  """
  trimmomatic ${params.read_type} \
              ${fq_reads[0]} \
              ${fq_reads[1]} \
              ${fq_reads[0]}_paired \
              ${fq_reads[0]}_unpaired \
              ${fq_reads[1]}_paired \
              ${fq_reads[1]}_unpaired \
              ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:True \
              LEADING:3 TRAILING:3 MINLEN:36
  """
}
