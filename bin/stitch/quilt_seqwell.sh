#!/bin/bash
#SBATCH -p compute
#SBATCH -q batch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 30
#SBATCH --mem 256G
#SBATCH --time 0-24:00
################################################################################
# Try to get QUILT running with the "good" seqWell data. This is a test to see
# if I can get the file formats correct. Only use Chr 1 reference data.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-23
################################################################################

##### VARIABLES #####

BASE_DIR=/projects/compsci/vmp/lcgbs_ssif

SCRATCH_DIR=/fastscratch/widmas

BAM_SRC_DIR=/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_seqwell_NovaSeq_full/bams

BAM_DEST_DIR=${SCRATCH_DIR}/bams

QUILT=${BASE_DIR}/singularity/quilt_1.0.4.sif

##### MAIN #####

# Load Singularity module.
module load singularity

# Go to scratch directory.
cd ${SCRATCH_DIR}

# Copy BAMs to fastscratch.
mkdir -p ${BAM_DEST_DIR}
cp ${BAM_SRC_DIR}/*.bam ${BAM_DEST_DIR}
cp ${BAM_SRC_DIR}/*.bai ${BAM_DEST_DIR}

# Make a file that lists the BAM filenames.
ls -1d ${BAM_DEST_DIR}/*.bam > bamfiles.txt

# Run QUILT.
singularity run ${QUILT} /opt/QUILT/QUILT.R \
   --chr=1 \
   --outputdir=quilt_seqwell \
   --bamlist=bamfiles.txt \
   --reference_haplotype_file=sanger_chr1_do_snps.hap.gz \
   --reference_legend_file=sanger_chr1_do_snps.legend.gz \
   --reference_sample_file=sanger_chr1_do_snps.samples \
   --nGen=40 \
   --nCores=20 \
   --addOptimalHapsToVCF=TRUE \
   --record_interim_dosages=TRUE \
   --save_prepared_reference=TRUE


#   --genetic_map_file=chr1_gen_map.txt \ # This takes hours to run as QUILT validates the map.
#   --make_plots=TRUE \ # This caused QUILT to crash. But I can make plots in R using the QUILT container.
#   --plot_per_sample_likelihoods=TRUE \

cp ${SCRATCH_DIR}/quilt_seqwell/quilt.1.vcf.gz /projects/compsci/vmp/lcgbs_ssif/results/quilt/quilt_do_seqwell.1.vcf.gz
cp ${SCRATCH_DIR}/quilt_seqwell/quilt.1.vcf.gz.tbi /projects/compsci/vmp/lcgbs_ssif/results/quilt/quilt_do_seqwell.1.vcf.gz.tbi

