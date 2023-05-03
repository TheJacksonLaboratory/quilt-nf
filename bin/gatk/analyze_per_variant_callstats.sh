#!/bin/bash
#SBATCH --job-name=pv_stats_plots
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 12:00:00
#SBATCH --mem=250G
#SBATCH --ntasks=1

# take gvcf directory from pipeline run as only argument
# i.e. data/4WC_ddRADseq_NovaSeq/gvcfs/

containerDir=/projects/compsci/vmp/USERS/widmas/stitch-nf/singularity_cache
singularity exec ${containerDir}/lcgbs_hr_qtl2_et_al.sif Rscript bin/gatk/per_variant_callstats.R $1
