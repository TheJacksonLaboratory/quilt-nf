process DOWNSAMPLE_BAM {

  cpus 1
  memory 50.GB
  time '04:00:00'
  //errorStrategy 'retry' 

  container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_post_downsample_coverage.txt", mode:'copy'

  input:
  tuple val(sampleID), val(coverage), file(bam)

  output:
  tuple file("*_downsampled.bam"), emit: downsampled_bam
  tuple file("*_post_downsample_coverage.txt"), emit: downsampled_depth_out
  
  script:
  log.info "----- Downsampling Reads from Sample: ${sampleID} -----"

  """
  echo ${coverage}
  echo ${params.downsample_to_cov}
  downsampling_coef=\$(awk 'BEGIN { print ( ${params.downsample_to_cov} / ${coverage} ) }')
  echo \$downsampling_coef

  function check_downsampling() {

  if awk -v num="\$1" 'BEGIN{ if (num > 1) exit 0; else exit 1; }'; then
    echo "Downsampling coefficient is greater than 1. Taking all reads."  

    # send the original bam
    samtools view -b ${bam} > ${sampleID}_skipped_downsampled.bam
    samtools index ${sampleID}_skipped_downsampled.bam

    # calculate new genome-wide coverage to verify
    samtools depth ${sampleID}_skipped_downsampled.bam -a | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_post_downsample_coverage.txt

  else
    echo "Downsampling coefficient is less than 1. Downsampling the .bam file."
    
    # downsample the bam file
    samtools view -b -s \$downsampling_coef ${bam} > ${sampleID}_downsampled.bam
    samtools index ${sampleID}_downsampled.bam

    # calculate new genome-wide coverage to verify
    samtools depth ${sampleID}_downsampled.bam -a | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_post_downsample_coverage.txt

    
  fi
  }

  check_downsampling "\$downsampling_coef"
  """


}
