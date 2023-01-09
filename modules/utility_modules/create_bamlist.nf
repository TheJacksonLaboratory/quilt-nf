process CREATE_BAMLIST {

  cpus 1
  memory 15.GB
  time '00:30:00'

  // publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/concatenated_reads' : 'concatenated_reads' }", pattern: "*fastq.gz", mode:'copy'

  input:
  tuple val(bams)

  output:
  tuple val(bamlist), emit: bam_list

  script:
  log.info "----- Create List of .bam Files for STITCH -----"

  """
  cat ${bams} > STITCH_bamlist.txt
  """
}