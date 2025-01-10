#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=MAKE_QUILT_REFERENCE-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.out

cd $SLURM_SUBMIT_DIR

LCGBS_DIR=/projects/compsci/vmp/lcgbs_ssif

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/23.10.1

# RUN PIPELINE
nextflow main.nf \
--workflow make_quilt_reference_data \
-profile sumner2 \
-w '/flashscratch/widmas/MAKE_QUILT_REFERENCE_DATA/work' \
--pubdir ${LCGBS_DIR}/results/make_quilt_reference_data \
--cross_type 'do' \
-resume