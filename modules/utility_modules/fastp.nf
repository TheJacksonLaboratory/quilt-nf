process FASTP {

  tag "$sampleID"

  cpus 1
  memory 50.GB
  time '3:00:00'

  container 'docker://sjwidmay/fastp_nf:fastp'
  
  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*filtered_trimmed*"), file("*_fastp_report.html"), emit: fastp_filtered

  script:

  if (params.library_type == "ddRADseq"){
   dedup_mode="--dont_eval_duplication"
  }
  if (params.library_type == "seqwell"){
   dedup_mode="-D"
  }

  """
  /fastp -i ${fq_reads[0]} \\
	-I ${fq_reads[1]} \\
	-o ${sampleID}_R1_filtered_trimmed.fq \\
	-O ${sampleID}_R2_filtered_trimmed.fq \\
	--detect_adapter_for_pe \\
	-g \\
	-c \\
	$dedup_mode \\
	-p \\
	--html ${sampleID}_fastp_report.html
  """
}
