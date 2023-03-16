process PILEUPS_TO_BAM {

  cpus 1
  memory 30.GB
  time '1:00:00'

  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

  input:
  tuple val(sampleID), file(bam), file(mpileup), file(coverage)

  output:
  tuple val(sampleID), file("*_covered.bam"), emit: filtered_bam

  script:
  log.info "----- Filtering Alignment from ${sampleID} to Pileup Sites -----"

  """
  awk '\$4 > 5 {print \$1"\t"\$2"\t"\$2}' ${mpileup} > ${sampleID}.bed
  bedtools intersect -abam ${bam} -b ${sampleID}.bed > ${sampleID}_covered.bam
  """
}