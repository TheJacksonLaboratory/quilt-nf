
process FASTP {

  tag "$sampleID"

  cpus 1
  memory 50.GB
  time '3:00:00'

  container 'docker://sjwidmay/fastp_nf:fastp'
  // container 'quay.io/biocontainers/fastp:0.23.2--h5f740d0_3'

  //publishDir "${params.sample_folder}/fastp", pattern:"*_fastp_report.html", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*filtered_trimmed*"), file("*_fastp_report.html"), emit: fastp_filtered

  script:
  log.info "----- FASTP Running on: ${sampleID} -----"

  //if (params.read_type == "SE"){
  //  mode_HQ="-S -M"
  //  inputfq="${fq_reads[0]}"
  //}
  //if (params.read_type == "PE"){
  //  mode_HQ="-M"
  //  inputfq="${fq_reads[0]} ${fq_reads[1]}"
  //}

  """
  /fastp -i ${fq_reads[0]} \\
	-I ${fq_reads[1]} \\
	-o ${sampleID}_R1_filtered_trimmed.fq \\
	-O ${sampleID}_R2_filtered_trimmed.fq \\
	--detect_adapter_for_pe \\
	-g \\
	-c \\
	-D \\
	-p \\
	--html ${sampleID}_fastp_report.html
  """
}
