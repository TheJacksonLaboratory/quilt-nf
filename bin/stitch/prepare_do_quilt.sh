#!/bin/bash
#SBATCH --qos dev
#SBATCH --partition dev
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 32G
#SBATCH --time 0-1:00
################################################################################
# Prepare data files for QUILT. Using DO seqWell samples from the "good" run.
# Use Chr 1 only.
#
# Make a reference data file first.
# This will involve preparing the following files:
# 1. Reference haplotype file in IMPUTE format.
# 2. Reference legend file in IMPUTE format.
# 3. Genetic and physical map.
#
# The QUILT docs say: "Note that QUILT requires a reference panel in the IMPUTE 
# .hap and .legend format. One easy way to make these from a haplotype VCF is by 
# using the bcftools convert --haplegendsample command."
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-23
################################################################################

##### VARIABLES #####

# Scripts directory.
SCRIPT_DIR='/projects/compsci/vmp/lcgbs_ssif/scripts'

# Containers.
BCFTOOLS=~/containers/samtools_1.10.sif

R=~/containers/bioconductor.sif

QUILT=~/projects/compsci/vmp/lcgbs_ssif/singularity/quilt_1.0.4.sif

##### MAIN #####

module load singularity

# Filter the Sanger VCF to retain SNPs on Chr 1.
bash filter_sanger_vcf.sh

# Create genetic map file.
singularity exec ${R} Rscript ${SCRIPT_DIR} make_quilt_map_file.R

# Filter the Sanger VCF.
bash ${SCRIPT_DIR}/filter_sanger_vcf.sh



