#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=STITCH-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN PIPELINE
nextflow main.nf \
--workflow stitch \
-profile sumner \
--sample_folder '/projects/compsci/vmp/USERS/widmas/stitch-nf/data/CC_data' \
--gen_org mouse \
--do_mice \
--pubdir '/fastscratch/STITCH_outputDir' \
-w '/fastscratch/STITCH_outputDir/work' \
--nFounders 8 \
--comment "This script will run haplotype inference using STITCH on DO plexWell data"
