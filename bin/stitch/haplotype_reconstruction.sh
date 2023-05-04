#!/bin/bash
#SBATCH --nodes 1 # number of nodes
#SBATCH --ntasks 30 # number of cores
#SBATCH --mem 256G # memory pool for all cores
#SBATCH -t 0-8:00 # time (D-HH:MM)

GENOME=grcm39

CONTAINER=~/containers/r_qtl2.sif

RSCRIPT=/projects/compsci/vmp/lcgbs_ssif/scripts/haplotype_reconstruction.R

module load singularity

singularity exec ${CONTAINER} Rscript ${RSCRIPT} ${GENOME}


# --qos dev
# --partition dev
