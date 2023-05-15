#!/usr/bin/env Rscript

# Load required packages
library(QUILT)

# Set arguments
args <- commandArgs(trailingOnly = TRUE)

#args <- c("/fastscratch/STITCH_outputDir/work/97/a0c11849aadfe4a0463064065e590a/bamlist.txt",
#          7,
#          "/fastscratch/STITCH_outputDir/work/48/260f7a481bbd78679dbbb46040613c/sanger_chr7_do.hap.gz",
#          "/fastscratch/STITCH_outputDir/work/48/260f7a481bbd78679dbbb46040613c/sanger_chr7_do.samples",
#          "/fastscratch/STITCH_outputDir/work/48/260f7a481bbd78679dbbb46040613c/sanger_chr7_do.legend.gz")

#args[1] = bamlist (from channel)
#args[2] = chromosome number (from channel)
#args[3] = reference haplotype file (if DO)
#args[4] = reference sample file (if DO)
#args[5] = reference legend file (if DO)

# Input files
# Path to original bamlist
mouse_bamlist <- args[1]

# Chromosome to be analyzed
mouse_chr <- paste0(args[2])

QUILT::QUILT(chr = mouse_chr,
               regionStart = 5000000,
               regionEnd = 10000000,
               buffer = 10000,
               bamlist = mouse_bamlist,
               outputdir = paste0(getwd(), "/"),
               reference_haplotype_file = args[3],
               reference_sample_file = args[4],
               reference_legend_file = args[5],
               nGen = 41,
               nCores = 20,
               addOptimalHapsToVCF=TRUE,
               record_interim_dosages=TRUE,
               save_prepared_reference=TRUE)
