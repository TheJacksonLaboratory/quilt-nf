#!/bin/bash

#SBATCH --job-name=subsample_CC
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 10:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1

homeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf
fastqDir=${homeDir}/data/CC_data


# Fire up the singularity container hopefully
singularity pull docker://quay.io/biocontainers/seqtk:1.3--hed695b0_2

# identify read pairs with one set of files
files=${fastqDir}/GES15*R1*
for f in ${files}
do
    # find the sample name
    sample=$(echo ${f} | cut -f 1 -d "." | cut -f 10 -d "/" | cut -f 1-5 -d "_")
    echo ${sample}
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/${sample}_R1.fastq.gz > ${fastqDir}/${sample}_R1.fa
    singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -a ${fastqDir}/${sample}_R2.fastq.gz > ${fastqDir}/${sample}_R2.fa

    # subsample CC reads
    for i in {1..5}
        do
        # generate seed - this is super important to distinguish actual subsamples
        seed=$(echo $(( $RANDOM % 100 + 1 )))

        # sample from fasta
        singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s${seed} ${fastqDir}/${sample}_R1.fa 2000000 > ${fastqDir}/sub${i}_${sample}_R1.fa
        singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk sample -s${seed} ${fastqDir}/${sample}_R2.fa 2000000 > ${fastqDir}/sub${i}_${sample}_R2.fa

        # convert back to fastq
        singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F 'F' ${fastqDir}/sub${i}_${sample}_R1.fa > ${fastqDir}/sub${i}_${sample}_R1.fastq
        singularity exec ${homeDir}/seqtk_1.3--hed695b0_2.sif seqtk seq -F 'F' ${fastqDir}/sub${i}_${sample}_R2.fa > ${fastqDir}/sub${i}_${sample}_R2.fastq
        
        # compress
        gzip ${fastqDir}/sub${i}_${sample}_R1.fastq --force
        gzip ${fastqDir}/sub${i}_${sample}_R2.fastq --force
        chmod 775 ${fastqDir}/sub${i}_${sample}_R1.fastq.gz
        chmod 775 ${fastqDir}/sub${i}_${sample}_R2.fastq.gz
        done
done
