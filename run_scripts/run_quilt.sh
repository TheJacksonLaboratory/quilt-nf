#!/bin/bash
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --job-name=QUILT-NF-LCWGS
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
--sample_folder ${PATH/TO/FASTQS} \
--pubdir ${PATH/TO/RESULTS_DIRECTORY} \
--extension 'fastq.gz' \
--pattern="*_R{1,2}*" \
--library_type "seqwell" \
--run_name $1 \ # results will appear in pubdir/run_name
-w ${PATH/TO/PROCESS_WORKING_DIRECTORIES} \
--bin_shuffling_file ${PATH/TO/BIN_SHUFFLE_FILE} \
--downsample_to_cov ${PATH/TO/DOWNSAMPLING_VALUE_FILE} \
--cross_type 'do' \
--align_only \
--covar_file${PATH/TO/QTL2_COVARFILE} \
--comment "This script will run haplotype inference on lcWGS data from complex mouse crosses" \
-resume
