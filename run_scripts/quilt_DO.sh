#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 36:00:00
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
--sample_folder '/projects/compsci/vmp/lcgbs_ssif/data/novaseq_seqwell_full/DO' \
--gen_org mouse \
--pubdir '/projects/compsci/vmp/lcgbs_ssif/results/quilt' \
--extension 'fastq.gz' \
--pattern="*_R{1,2}*" \
--library_type "seqwell" \
--run_name $1 \
-w '/flashscratch/STITCH_outputDir/work' \
--downsample_to_cov '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/downsampling_values_small.csv' \
--cross_type 'do' \
--ref_file_dir '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders' \
--covar_file '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/DO_covar.csv' \
--comment "This script will run haplotype inference on DO lcWGS data" \
-resume
