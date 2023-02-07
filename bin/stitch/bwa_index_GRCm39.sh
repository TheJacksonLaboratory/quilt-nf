#!/bin/bash
#SBATCH -J BWA_INDEX_GRCm39
#SBATCH --mem 30GB
#SBATCH -t 2:30:00
#SBATCH -p compute

singularity exec /projects/omics_share/meta/containers/quay.io/biocontainers/bwakit:0.7.17.dev1--hdfd78af_1 \ 
        bwa index -a bwtsw \ 
        /projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.primary_assembly.fa \ 
        Mus_musculus.GRCm39.dna.toplevel.fa