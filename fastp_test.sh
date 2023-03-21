#!/bin/bash

#SBATCH --job-name=fastp_test
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 01:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1

read1=/fastscratch/STITCH_outputDir/work/e7/1ffbfb7263b926b9a0c400f77ddab6/plexWell-M24_GT23-00380_CTGGTCGT-ATCTATTT_S36_L001_R1_001.fastq.gz_filtered_trimmed
read2=/fastscratch/STITCH_outputDir/work/e7/1ffbfb7263b926b9a0c400f77ddab6/plexWell-M24_GT23-00380_CTGGTCGT-ATCTATTT_S36_L001_R2_001.fastq.gz_filtered_trimmed
outDir=/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_seqwell_NovaSeq/fastp_test
singularity exec docker://sjwidmay/fastp_nf:fastp /fastp -i ${read1} -I ${read2} -o ${outDir}/test.R1.fq -O ${outDir}/test.R2.fq --detect_adapter_for_pe -g -D -c -p

# polyG trimming
#        -g \   
