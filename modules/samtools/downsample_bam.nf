process DOWNSAMPLE_BAM {

  cpus 1
  memory 50.GB
  time '01:00:00'
  //errorStrategy 'retry' 
  maxRetries 3

  container 'quay.io/biocontainers/samtools:1.16.1--h00cdaf9_2'

  publishDir "${params.pubdir}/${params.run_name}/coverage", pattern:"*_downsampled_coverage.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(coverage_file)

  output:
  tuple val(sampleID), file("*_downsampled.bam"), emit: downsampled_bam
  tuple file("*_downsampled_coverage.txt"), emit: downsampled_depth_out

  script:
  log.info "----- Downsampling Reads from Sample: ${sampleID} -----"

  """
  # pull genome-wide coverage from previous process
  coverage=$(cat ${coverage_file})
  coverage=$(echo "$coverage" | bc -l)
  
  # desired final coverage of bam file
  final_coverage=${params.downsampleToCov}

  # determine final downsampling coefficient
  downsampling_coef=$(echo "scale=5; $final_coverage / $coverage" | bc)
  
  # if the downsampling coefficient is greater than 1, just take the original bam
  if (( $(echo "$downsampling_coef > 1" | bc -l) )); the
    echo "Downsampling coefficient is greater than 1. Taking all reads."
    
    # rename bam and send with the others
    mv ${bam} > ${sampleID}_skipped_downsampled.bam

    # send original coverage along
    mv ${coverage_file} > ${sampleID}_skipped_downsampled_coverage.txt
  
  else
    echo "Downsampling coefficient is less than 1. Downsampling the .bam file."

    # downsample the bam file
    samtools view -b -s ${downsampling_coef} ${bam} > ${sampleID}_downsampled.bam

    # calculate new genome-wide coverage to verify
    samtools depth ${sampleID}_downsampled.bam -a | awk 'BEGIN{sum=0} {sum += \$3; n++} END{print sum/n}' > ${sampleID}_downsampled_coverage.txt
  fi
  """
}
