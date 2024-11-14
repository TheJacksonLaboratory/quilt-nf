process QUILT_LD_PRUNING {
  tag "$chr"

  cpus 1
  memory 50.GB
  time '01:00:00'

  container 'docker://sjwidmay/quilt_nf:plink_bcftools'

  publishDir "${params.pubdir}/${params.run_name}/quilt_vcfs", pattern:"*", 
mode:'copy'

  input:
  tuple val(chr), file(sample_genos), file(sample_genos_index)

  output:
  tuple val(chr), file("pruned.quilt.*.vcf.gz"), file("pruned.quilt.*.vcf.gz.csi"), emit: pruned_quilt_vcf

  script:
  log.info "----- Performing LD Pruning on QUILT Variants for Chromosome: ${chr} -----"

  """
  # find markers in LD
  plink --vcf ${sample_genos} \
    --double-id \
    --maf 0.01 \
    --geno \
    --set-missing-var-ids @:# \
    --chr ${chr} \
    --indep-pairwise 50 20 0.9 \
    --make-bed

  # get a tab separated marker list for bcftools to read
  awk '{sub(/:/, "\\t"); print}' plink.prune.in  > marker.list.txt

  # perform LD pruning
  bcftools view ${sample_genos} \
    -R marker.list.txt \
    --output-type z \
    -o pruned.quilt.${chr}.vcf.gz
  
  # index new vcf
  bcftools index pruned.quilt.${chr}.vcf.gz
  """
}
