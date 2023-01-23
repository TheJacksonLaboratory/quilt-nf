#!/bin/bash

#SBATCH --job-name=subsample_CC
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 10:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1


fastqDir=/fastscratch/widmas/CC_fastqs
homeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf

# Fire up the singularity container hopefully
singularity pull docker://quay.io/biocontainers/seqtk:1.3--hed695b0_2

files=${fastqDir}/GES15-067*
for f in ${files}
do
    # find the sample name
    # echo $f |  cut -f 1 -d "." | cut -f 5 -d "/"
    sample=$(echo ${f} | cut -f 1 -d "." | cut -f 5 -d "/")
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${f} > ${fastqDir}/${sample}.fa

    # subsample CC reads
    for i in {1..20}
        do
        singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample ${fastqDir}/${sample}.fa 2000000 > ${homeDir}/test/wgs/mouse/sub${i}_${sample}.fa
        singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F '#' ${homeDir}/test/wgs/mouse/sub${i}_${sample}.fa > ${homeDir}/test/wgs/mouse/sub${i}_${sample}.fastq
        gzip ${homeDir}/test/wgs/mouse/sub${i}_${sample}.fastq --force
        done
done
