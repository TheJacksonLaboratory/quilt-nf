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
singularity pull docker://sjwidmay/stitch_nf:bbtools

# identify read pairs with one set of files
files=${fastqDir}/GES15*R1*
for f in ${files}
do
    # find the sample name
    sample=$(echo ${f} | cut -f 1 -d "." | cut -f 10 -d "/" | cut -f 1-5 -d "_")
    echo ${sample}

    # subsample CC reads
    for i in {1..5}
        do
        echo ${i}
        # generate seed - this is super important to distinguish actual subsamples
        # seed=$(echo $(( $RANDOM % 100 + 1 )))

        # sample from fasta
        singularity exec ${homeDir}/stitch_nf_bbtools.sif reformat.sh in1=${fastqDir}/${sample}_R1.fastq.gz in2=${fastqDir}/${sample}_R2.fastq.gz out1=${fastqDir}/sub${i}_${sample}_R1.fastq.gz out2=${fastqDir}/sub${i}_${sample}_R2.fastq.gz samplereadstarget=5000000

        # convert back to fastq
        chmod 775 ${fastqDir}/sub${i}_${sample}_R1.fastq.gz
        chmod 775 ${fastqDir}/sub${i}_${sample}_R2.fastq.gz
        done
done
