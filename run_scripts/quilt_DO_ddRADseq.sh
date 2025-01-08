#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF-DDRAD
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.out

cd $SLURM_SUBMIT_DIR

LCGBS_DIR=/projects/compsci/vmp/lcgbs_ssif

QUILT_DIR=/projects/compsci/vmp/USERS/widmas/quilt-nf

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/23.10.1

# RUN PIPELINE
nextflow main.nf \
--workflow quilt \
-profile sumner2 \
--sample_folder ${LCGBS_DIR}/data/novaseq_ddradseq2/stacks_fastqs_w_umis/DO \
--gen_org mouse \
--pubdir ${LCGBS_DIR}/results/quilt \
--extension 'fq.gz' \
--pattern "*.{1,2}.*" \
--library_type "ddRADseq" \
--run_name $1 \
-w '/flashscratch/widmas/QUILT/work' \
--bin_shuffling_file ${QUILT_DIR}/data/shuffle_bins.csv \
--downsample_to_cov ${QUILT_DIR}/data/downsampling_values.csv \
--cross_type 'do' \
--align_only \
--covar_file ${QUILT_DIR}/data/DO_covar.csv \
--comment "This script will run haplotype inference on DO ddRADseq data" \
-resume
