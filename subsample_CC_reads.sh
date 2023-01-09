#!/bin/bash

#SBATCH --job-name=subsample_CC
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 03:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1


fastqDir=/fastscratch/widmas/CC_fastqs
homeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf
CC007=GES15-06702-CC007-Unc_GAGTGGAT_HMGWYCCXX_L007_001
CC062=GES15-06751-CC062-Unc_ATGTCAGA_HMHGCCCXX_L008_001

# Fire up the singularity container hopefully
singularity pull docker://quay.io/biocontainers/seqtk:1.3--hed695b0_2

singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/${CC007}_R1.fastq.gz > ${fastqDir}/${CC007}_R1.fa
singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/${CC007}_R2.fastq.gz > ${fastqDir}/${CC007}_R2.fa
singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/${CC062}_R1.fastq.gz > ${fastqDir}/${CC062}_R1.fa
singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/${CC062}_R2.fastq.gz > ${fastqDir}/${CC062}_R2.fa

# loop through 20 times and sample 50k reads
for i in {1..20}
do
    # subsample CC007 reads
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s20 ${fastqDir}/${CC007}_R1.fa 50000 > ${homeDir}/test/wgs/mouse/CC007-Unc_50k_${i}_R1.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F '#' ${homeDir}/test/wgs/mouse/CC007-Unc_50k_${i}_R1.fa > ${homeDir}/test/wgs/mouse/${CC007}_subs${i}_R1.fastq 
    gzip ${homeDir}/test/wgs/mouse/${CC007}_subs${i}_R1.fastq

    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s20 ${fastqDir}/${CC007}_R2.fa 50000 > ${homeDir}/test/wgs/mouse/CC007-Unc_50k_${i}_R2.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F '#' ${homeDir}/test/wgs/mouse/CC007-Unc_50k_${i}_R2.fa > ${homeDir}/test/wgs/mouse/${CC007}_subs${i}_R2.fastq 
    gzip ${homeDir}/test/wgs/mouse/${CC007}_subs${i}_R2.fastq

    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s20 ${fastqDir}/${CC062}_R1.fa 50000 > ${homeDir}/test/wgs/mouse/CC062-Unc_50k_${i}_R1.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F '#' ${homeDir}/test/wgs/mouse/CC062-Unc_50k_${i}_R1.fa > ${homeDir}/test/wgs/mouse/${CC062}_subs${i}_R1.fastq 
    gzip ${homeDir}/test/wgs/mouse/${CC062}_subs${i}_R1.fastq

    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s20 ${fastqDir}/${CC062}_R2.fa 50000 > ${homeDir}/test/wgs/mouse/CC062-Unc_50k_${i}_R2.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F '#' ${homeDir}/test/wgs/mouse/CC062-Unc_50k_${i}_R2.fa > ${homeDir}/test/wgs/mouse/${CC062}_subs${i}_R2.fastq 
    gzip ${homeDir}/test/wgs/mouse/${CC062}_subs${i}_R2.fastq
done


