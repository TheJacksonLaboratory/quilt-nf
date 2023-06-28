#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q long
#SBATCH -t 96:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN PIPELINE
nextflow main.nf \
--workflow quilt \
-profile sumner \
--sample_folder '/projects/compsci/vmp/USERS/widmas/stitch-nf/data/4WC_seqwell_NovaSeq_full' \
--gen_org mouse \
--pubdir '/projects/compsci/vmp/lcgbs_ssif/results/quilt' \
--run_name $1 \
-w '/fastscratch/STITCH_outputDir/work' \
--downsample_to_cov 30 \
--cross_type 'genail4' \
--ref_file_dir '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/' \
--covar_file '/projects/compsci/vmp/lcgbs_ssif/data/GigaMUGA/4WC_covar.csv' \
--comment "This script will run haplotype inference on 4WC lcWGS data" \
-resume
