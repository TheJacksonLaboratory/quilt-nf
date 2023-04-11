process COMBINE_GVCF {
// Note about this module
  tag "$sampleID"

  cpus = 4
  memory {100.GB * task.attempt}
  errorStrategy 'retry' 
  maxRetries 4
  time = '05:30:00'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(chrom), val(sampleID), file(gvcf)

  output:
  tuple val(chrom), file("*_combined.g.vcf.gz"), emit: chr_vcf
  //tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- Combining Chromosome ${chrom} GVCFs -----"

  """
  ls ${gvcf} | sed 's/^/-V /' > gvcf.list
  sed '1,'"$(( $(wc -l < gvcf.list) - 1 ))"' s/$/ \\/g' gvcf.list > gvcf.list.final
  cat ${projectDir}/bin/gatk/combineGVCF_template.sh && echo && cat gvcf.list.final > chr${chrom}_combine_gvcfs.sh
  bash chr${chrom}_combine_gvcfs.sh ${params.ref_fa} chr${chrom}_combined.g.vcf.gz
  """
}