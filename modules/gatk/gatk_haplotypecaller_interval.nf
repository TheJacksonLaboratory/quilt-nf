process GATK_HAPLOTYPECALLER_INTERVAL {
// Note about this module
  tag "$sampleID"

  cpus = 4
  memory {100.GB * task.attempt}
  errorStrategy 'retry' 
  maxRetries 4
  time = '05:30:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.sample_folder}/gvcfs", pattern: "*.g.vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), file(bam), file(bai), val(chrom)

  output:
  tuple val(chrom), val(sampleID), file("*_genotyped.vcf.gz"), emit: vcf
  // tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- GATK Haplotype Caller Running on Chromosome ${chrom} for sample: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_HaplotypeCaller_${chrom}.g.vcf.gz \
  -L ${chrom} \
  -ERC GVCF \
  -stand-call-conf 30 \
  --max-num-haplotypes-in-population ${params.nFounders}

  gatk --java-options "-Xmx${my_mem}G" GenotypeGVCFs  \
  -R ${params.ref_fa} \
  -V ${sampleID}_HaplotypeCaller_${chrom}.g.vcf.gz \
  -O ${sampleID}_chr${chrom}_genotyped.vcf.gz \
  """
}
