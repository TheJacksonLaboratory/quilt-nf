#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow/23.10.1

# RUN PIPELINE
nextflow main.nf \
--workflow quilt \
-profile sumner2 \
--sample_folder '/projects/compsci/vmp/lcgbs_ssif/data/novaseq_ddradseq2/stacks_fastqs/raw/DO' \
--gen_org mouse \
--pubdir '/projects/compsci/vmp/lcgbs_ssif/results/quilt' \
--extension 'fq.gz' \
--pattern "*.{1,2}.*" \
--library_type "ddRADseq" \
--run_name $1 \
-w '/flashscratch/STITCH_outputDir/work' \
--bin_shuffling_file '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/shuffle_bins.csv' \
--downsample_to_cov '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/downsampling_values.csv' \
--cross_type 'do' \
--ref_file_dir '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders' \
--covar_file '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/DO_covar.csv' \
--comment "This script will run haplotype inference on DO ddRADseq data" \
-resume
