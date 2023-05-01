#!/bin/bash
#SBATCH --job-name=lcgbs_pv_stats
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 2:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1

homeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf
gvcfDir=${homeDir}/data/4WC_ddRADseq_NovaSeq/gvcfs
containerDir=${homeDir}/singularity_cache

chr_gvcfs=${gvcfDir}/chr*
#gvcf_test=${gvcfDir}/chr1_genotyped.vcf.gz

for g in ${chr_gvcfs}
do
    chrom=$(echo ${g} | cut -f 1 -d "." | cut -f 11 -d "/")
    
    echo ${chrom}
    singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view -m2 -M2 -v snps ${g} | \
    singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query --print-header -f '%CHROM\t%POS\t%REF\t%ALT\t%DP[\t%AD]\n' | \
    sed 's/[[# 0-9]*]//g' > ${gvcfDir}/${chrom}_per_variant_depth_stats.txt
done
