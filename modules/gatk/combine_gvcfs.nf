process COMBINE_GVCF {
// Note about this module
  tag "$chrom"

  cpus = 4
  memory 101.GB
  time = '01:30:00'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(chrom), val(sampleID), file(gvcf)

  output:
  tuple val(chrom), file("*.txt"), emit: chr_vcf
  //tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- Combining Chromosome ${chrom} GVCFs -----"

  """
  
  ls ${gvcf} > gvcf.list.txt 
  
  sed 's/^/-V /' gvcf.list > gvcf.list.2
  sed '1,'"\$(( \$(wc -l < gvcf.list.2) - 1 ))" 's/\$/ \\\/g' gvcf.list.2 > gvcf.list.final
  cat ${projectDir}/bin/gatk/combineGVCF_template.sh && echo && cat gvcf.list.final > chr${chrom}_combine_gvcfs.sh
  bash chr${chrom}_combine_gvcfs.sh ${params.ref_fa} chr${chrom}_combined.g.vcf.gz

  """
}
