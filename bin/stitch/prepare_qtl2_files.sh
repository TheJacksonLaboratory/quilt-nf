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

BASE_DIR=/projects/compsci/vmp/USERS/widmas/stitch-nf

# Founder VCF.
FOUNDER_FILE=/fastscratch/widmas/sanger_chr1_do_snps.vcf.gz

# Sample VCF (from QUILT).
SAMPLE_VCF=/fastscratch/widmas/quilt_seqwell/quilt.1.vcf.gz

# Sample metadat file.
META_FILE=/projects/compsci/vmp/USERS/widmas/lcGBS_wf/data/DO_covar.csv

# SNP annotation file.
SNP_FILE=/fastscratch/widmas/chr1_gen_map.txt

# qtl2 output directory
QTL2_DIR=/fastscratch/widmas/qtl2_do_seqwell

# R script to create files.
R_SCRIPT=${BASE_DIR}/bin/stitch/prepare_do_qtl2_files.R

# R container for R script
R=/projects/compsci/vmp/lcgbs_ssif/singularity/bioconductor.sif


##### MAIN #####

module load singularity

# Create qtl2 output directory.
mkdir -p ${QTL2_DIR} 

# Run R script in container.
singularity exec ${R} Rscript ${R_SCRIPT} ${FOUNDER_FILE} ${SAMPLE_VCF} ${META_FILE} ${SNP_FILE} ${QTL2_DIR}

