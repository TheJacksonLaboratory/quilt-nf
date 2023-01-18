#!/bin/bash

#SBATCH --job-name=subsample_CC
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 03:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1


fastqDir=/fastscratch/widmas/CC_fastqs
homeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf

# Fire up the singularity container hopefully
singularity pull docker://quay.io/biocontainers/seqtk:1.3--hed695b0_2

files=${fastqDir}/GES15-067*
for f in ${fastqDir}/${files}
do
    # find the sample name
    sample=echo $f |  cut -f 1 -d "."

    # subsample CC reads
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/$f > ${fastqDir}/${sample}.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s20 ${fastqDir}/${sample}.fa 2000000 > ${homeDir}/test/wgs/mouse/${sample}_${i}.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F '#' ${homeDir}/test/wgs/mouse/${sample}_${i}.fa > ${homeDir}/test/wgs/mouse/${sample}_${i}.fastq
    gzip ${homeDir}/test/wgs/mouse/${sample}_${i}.fastq
done