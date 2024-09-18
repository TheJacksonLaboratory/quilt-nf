#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 36:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1

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
--sample_folder ${LCGBS_DIR}/data/novaseq_ddradseq2/stacks_fastqs/raw/4WC \
--gen_org mouse \
--pubdir ${LCGBS_DIR}/results/quilt \
--extension 'fq.gz' \
--pattern "*.{1,2}.*" \
--library_type "ddRADseq" \
--run_name $1 \
-w '/flashscratch/widmas/QUILT_outputDir/work' \
--bin_shuffling_file ${QUILT_DIR}/data/shuffle_bins.csv \
--downsample_to_cov ${QUILT_DIR}/data/downsampling_values_4WC.csv \
--cross_type 'genail4' \
--ref_file_dir ${LCGBS_DIR}/data/4wc_founders \
--covar_file ${LCGBS_DIR}/data/GigaMUGA/4WC_covar_quilt.csv \
--align_only FALSE \
--comment "This script will run haplotype inference on 4WC ddRADseq data"
