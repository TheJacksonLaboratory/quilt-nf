#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32G
#SBATCH --time 0-1:00
################################################################################
# Filter the Sanger VCF to include only the DO founders and high-quality, 
# biallelic SNPs on Chr 1.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-23
################################################################################


##### VARIABLES #####

BASE_DIR=/fastscratch/dgatti

# Full Sanger VCF with all samples. GRCm39.
SANGER_DIR=/projects/omics_share/mouse/GRCm39/genome/annotation/snps_indels/rel_2112_v8
SANGER_FILE=${SANGER_DIR}/mgp_REL2021_snps.vcf.gz

# GRCm39 reference FASTA file.
REF_FASTA=/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa

# Script directory for R scripts.
SCRIPT_DIR=/projects/compsci/vmp/lcgbs_ssif/scripts

# Output directory and filename.
OUT_DIR=/fastscratch/dgatti
OUT_FILE=${OUT_DIR}/sanger_hq_do_snps.vcf

# bcftools container
BCFTOOLS=~/containers/samtools_1.10.sif

# Bioconductor container
BIOCONDUCTOR=~/containers/bioconductor.sif


##### MAIN #####

# Filter the Sanger VCF to contain only Chr 1 SNPs with a PASS filter & write out
# to a temporary file.
singularity exec ${BCFTOOLS} bcftools filter \
       --regions 1 \
       --include 'INFO/INDEL=0 && FILTER="PASS" && TYPE="snp"' \
       --output-type z \
       --output ${OUT_DIR}/s_file1.vcf.gz \
       ${SANGER_FILE}

# Filter the Chr 1 VCF to include only SNPs from DO founders.
# Filter to include only polymorphic, biallelic SNPs.
singularity exec ${BCFTOOLS} bcftools view \
       --samples A_J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ \
       --genotype hom \
       --min-alleles 2 \
       --max-alleles 2 \
       --min-ac 2 \
       --output-type z \
       --output-file ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz \
       ${OUT_DIR}/s_file1.vcf.gz

# Index VCF.
singularity exec ${BCFTOOLS} bcftools index ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz

rm ${OUT_DIR}/s_file1.vcf.gz

# Filter the Sanger file further and create a C57BL/6J tab-delimited file.
# NOTE: THIS WRITE OUT A *.bgz VCF file. (VariantAnnotation default)
singularity exec ${BIOCONDUCTOR} Rscript ${SCRIPT_DIR}/filter_sanger_file.R ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz

# Remove the old VCF so that we don't get confused below. 
# We're working with the *.bgz file now.
rm ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz.csi

# Convert C57BL/6J tab-delimited file to VCF.
singularity exec ${BCFTOOLS} bcftools convert \
        --tsv2vcf ${OUT_DIR}/C57BL_6J.tab \
        --fasta-ref ${REF_FASTA} \
        --samples C57BL_6J \
        --output-type z \
        --output ${OUT_DIR}/C57BL_6J.vcf.gz

# Index the C57BL/6J file.
singularity exec ${BCFTOOLS} bcftools index ${OUT_DIR}/C57BL_6J.vcf.gz 

# Merge C57BL/6J into filtered Sanger VCF.
singularity exec ${BCFTOOLS} bcftools merge \
        --output-type z \
        --output ${OUT_DIR}/sanger_merged.vcf.gz \
        ${OUT_DIR}/C57BL_6J.vcf.gz \
        ${OUT_DIR}/sanger_chr1_do_snps.vcf.bgz

# Clean up and rename output VCF file.
# The founder VCF now has a *.gz extension. 
rm ${OUT_DIR}/C57BL_6J*
rm ${OUT_DIR}/sanger_chr1_do_snps.vcf.*

# Make the final file phased since all of the SNPs are homozygous.
zcat ${OUT_DIR}/sanger_merged.vcf.gz | sed '/^##/! s/\//\|/g' > ${OUT_DIR}/sanger_chr1_do_snps.vcf
singularity exec ${BCFTOOLS} bgzip ${OUT_DIR}/sanger_chr1_do_snps.vcf
rm ${OUT_DIR}/sanger_merged.vcf.gz

# Index the Sanger VCF.
singularity exec ${BCFTOOLS} bcftools index ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz

# Make the map file based on this VCF.
singularity exec ${BIOCONDUCTOR} Rscript ${SCRIPT_DIR}/make_quilt_map_file.R ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz 

# Convert VCF to IMPUTE2 format.
singularity exec ${BCFTOOLS} bcftools convert ${OUT_DIR}/sanger_chr1_do_snps.vcf.gz \
       --haplegendsample ${BASE_DIR}/sanger_chr1_do_snps

