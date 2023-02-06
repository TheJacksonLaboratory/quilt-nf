process RUN_STITCH {
  tag "$chr"
  
  memory 100.GB
  time '10:00:00'
  
  container 'docker://sjwidmay/stitch_nf:latest'
  container 'quay.io-biocontainers-bcftools-1.15--h0ea216a_2'

  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "stitch.*.txt", mode:'copy'
  publishDir "${params.sample_folder}/stitch_vcfs", pattern: "EM.all.*.RData", mode:'copy'
  
  input:
  tuple file(bamlist), val(chr), file(posfile)

  output:
  tuple val(chr), file("stitch.*.txt"), file("RData/EM.all.*.RData"), emit: stitch_output 

  script:
  log.info "----- Running STITCH on Chromosome ${chr} -----"

  """
  Rscript --vanilla ${projectDir}/bin/stitch/run_stitch.R ${bamlist} ${posfile} ${params.nFounders} ${chr}
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' stitch.${chr}.vcf.gz > stitch.${chr}.txt
  """
}
