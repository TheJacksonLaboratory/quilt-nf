#!/bin/bash

#SBATCH --job-name=subsample_CC
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 01:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1


fastqDir = /fastscratch/widmas/CC_fastqs/
homeDir = /projects/compsci/vmp/USERS/widmas/stitch-nf

# Fire up the singularity container hopefully
singularity run docker://quay.io/biocontainers/seqtk

# loop through 20 times and sample 50k reads
for i in 1:20 ;
do
    # subsample CC007 reads
     seqtk sample -s20 \ 
        ${fastqDir}/GES15-06702-CC007-Unc_GAGTGGAT_HMGWYCCXX_L007_001_R1.fastq.gz \ 
        50000 > ${homeDir}/test/wgs/mouse/CC007-Unc_50k_${i}_R1.fastq.gz
    
    seqtk sample -s20 \ 
        ${fastqDir}/GES15-06702-CC007-Unc_GAGTGGAT_HMGWYCCXX_L007_001_R2.fastq.gz \ 
        50000 > ${homeDir}/test/wgs/mouse/CC007-Unc_50k_${i}_R2.fastq.gz

    # subsample CC0062 reads
    seqtk sample -s20 \ 
        ${fastqDir}/GES15-06751-CC062-Unc_ATGTCAGA_HMHGCCCXX_L008_001_R1.fastq.gz \ 
        50000 > ${homeDir}/test/wgs/mouse/CC062-Unc_50k_${i}_R1.fastq.gz

    seqtk sample -s20 \ 
        ${fastqDir}/GES15-06751-CC062-Unc_ATGTCAGA_HMHGCCCXX_L008_001_R2.fastq.gz \ 
        50000 > ${homeDir}/test/wgs/mouse/CC062-Unc_50k_${i}_R2.fastq.gz
done


