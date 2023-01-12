#!/bin/bash
#SBATCH -J STITCH_TEST
#SBATCH --mem 30GB
#SBATCH -t 1:30:00
#SBATCH -p compute

singularity run docker://sjwidmay/stitch-nf:latest /projects/compsci/vmp/USERS/widmas/stitch-nf/bin/run_stitch.R /fastscratch/STITCH_outputDir/work/cc/73ed7fed40b46912004bb84482aa15/STITCH_bamlist.txt /fastscratch/STITCH_outputDir/work/cc/73ed7fed40b46912004bb84482aa15/STITCH_15_pos.txt 4 15
  