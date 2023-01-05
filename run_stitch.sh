#!/bin/bash

#SBATCH --job-name=stitch-nf-test
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN TEST PIPELINE
nextflow main.nf \
-profile sumner \
--workflow stitch \
--gen_org mouse \
--sample_folder '/fastscratch/widmas/CC_fastqs' \
--pubdir '/fastscratch/STITCH_outputDir' \
-w '/fastscratch/STITCH_outputDir/work'
