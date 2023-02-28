#!/bin/bash
singularity run /projects/compsci/omics_share/meta/containers/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img

bcftools view /projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8/mgp_REL2021_snps.vcf.gz \
-m2 -M2 -v snps \
-s C57BL_6NJ | \
sed 's/1\/1/0\/0/g' | \
sed 's/0\/1/0\/0/g' | \
sed 's/1\/0/0\/0/g' | \
sed 's/C57BL_6NJ/B/g' | bcftools view -Oz -o /fastscratch/widmas/B6.vcf.gz

# index B6 vcf
tabix -p vcf /fastscratch/widmas/B6.vcf.gz

# create the other main vcf
bcftools view /projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8/mgp_REL2021_snps.vcf.gz \
-m2 -M2 -v snps \
-s A_J,C57BL_6NJ,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ | \
sed 's/A_J/A/g' | \
sed 's/C57BL_6NJ/B/g' | \
sed 's/129S1_SvImJ/C/g' | \
sed 's/NOD_ShiLtJ/D/g' | \
sed 's/NZO_HlLtJ/E/g' | \
sed 's/CAST_EiJ/F/g' | \
sed 's/PWK_PhJ/G/g' | \
sed 's/WSB_EiJ/H/g' | \
bcftools view -Oz -o /fastscratch/widmas/notB6.vcf.gz

# index main vcf
tabix -p vcf /fastscratch/widmas/notB6.vcf.gz

# join the two, but select the altered B6 sample field
bcftools merge /fastscratch/widmas/B6.vcf.gz /fastscratch/widmas/notB6.vcf.gz --force-samples | \
bcftools view -s A,B,C,D,E,F,G,H -Oz -o bin/stitch/DO_founders.vcf.gz
