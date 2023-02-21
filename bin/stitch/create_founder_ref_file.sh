#!/bin/bash
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 12:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1

singularity exec /projects/compsci/omics_share/meta/containers/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view /projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8/mgp_REL2021_snps.vcf.gz --regions 19 -m2 -M2 -v snps -s A_J,C57BL_6NJ,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ | sed 's/1\/1/0\/0/2' | sed 's/0\/1/0\/0/2' | sed 's/1\/0/0\/0/2' | sed 's/A_J/A/g' | sed 's/C57BL_6NJ/B/g' | sed 's/129S1_SvImJ/C/g' | sed 's/NOD_ShiLtJ/D/g' | sed 's/NZO_HlLtJ/E/g' | sed 's/CAST_EiJ/F/g' | sed 's/PWK_PhJ/G/g' | sed 's/WSB_EiJ/H/g' > founders_GRCm39_letters.vcf.gz
