#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --output=%x.%j.out

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/23.10.1

# RUN PIPELINE
nextflow main.nf \
--workflow quilt \
-profile sumner2 \
--sample_folder '/projects/reinholdt-lab/DO_ESC/data/lcgbs' \
--gen_org mouse \
--pubdir '/projects/reinholdt-lab/DO_ESC/results/quilt' \
--extension 'fastq.gz' \
--pattern="*_R{1,2}.*" \
--library_type "seqwell" \
--run_name $1 \
-w '/flashscratch/widmas/QUILT/work' \
--bin_shuffling_file '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/DO_ESC_binshuffle.csv' \
--downsample_to_cov '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/downsampling_values.csv' \
--cross_type 'do' \
--ref_file_dir '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders' \
--covar_file '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/DO_ESC_covar.csv' \
--comment "This script will run haplotype inference on DO lcWGS data" \
--align_only \
-resume
