#!/bin/bash
#SBATCH --qos batch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 30
#SBATCH --mem 128G
#SBATCH --time 0-24:00
################################################################################
# Once QUILT has run, get the sample VCF and the founder VCF, subset to include
# only high-quality, biallelic, homozygous founder SNPs, and write out qtl2
# support files.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-23
################################################################################

##### VARIABLES #####

BASE_DIR=/projects/compsci/vmp/lcgbs_ssif

# Founder VCF.
FOUNDER_FILE=/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/4Founders.deepVar.vcf.gz

# Sample VCF (from QUILT).
SAMPLE_VCF=${BASE_DIR}/results/quilt/quilt_do_seqwell.1.vcf.gz

# Sample metadat file.
META_FILE=${BASE_DIR}/data/4WC_covar.csv

# SNP annotation file.
SNP_FILE=/fastscratch/dgatti/chr1_gen_map.txt

# qtl2 output directory
QTL2_DIR=/fastscratch/dgatti/qtl2_do_seqwell

# R script to create files.
R_SCRIPT=${BASE_DIR}/scripts/prepare_do_qtl2_files.R

# R container for R script
R=~/containers/bioconductor.sif


##### MAIN #####

module load singularity

# Create qtl2 output directory.
mkdir -p ${QTL2_DIR} 

# Run R script in container.
singularity exec ${R} R_script ${R_SCRIPT} ${FOUNDER_FILE} ${SAMPLE_FILE} ${META_FILE} ${SNP_FILE} ${QTL2_DIR}

