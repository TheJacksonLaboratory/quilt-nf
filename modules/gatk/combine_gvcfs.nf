process COMBINE_GVCF {
// Note about this module
  tag "$chrom"

  cpus = 4
  memory 100.GB
  time = '01:30:00'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(chrom), val(sampleID), file(gvcf), file(gvcf_index)

  output:
  tuple val(chrom), file("*combine_gvcfs.sh"), emit: chr_vcf
  //tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- Combining Chromosome ${chrom} GVCFs -----"

  """
  ls ${gvcf} > gvcf.list.txt
  sed 's/^/-V /' gvcf.list.txt > gvcf.list.V.txt
  nsamples_but_last=\$((\$(wc -l < gvcf.list.V.txt) - 1))
  sed "1,\${nsamples_but_last}s/\$/ \\\\\\/" gvcf.list.V.txt > gvcf.list.final.txt
  cat ${projectDir}/bin/gatk/combineGVCF_template.sh gvcf.list.final.txt > chr${chrom}_combine_gvcfs.sh
  bash chr${chrom}_combine_gvcfs.sh ${params.ref_fa} chr${chrom}_combined.g.vcf.gz
  """
}
