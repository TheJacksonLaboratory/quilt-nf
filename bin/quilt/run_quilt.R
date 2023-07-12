#!/usr/bin/env Rscript

# Load required packages
library(QUILT)

# Set arguments
args <- commandArgs(trailingOnly = TRUE)

#args[1] = bamlist (from channel)
#args[2] = chromosome number (from channel)
#args[3] = reference haplotype file (if DO)
#args[4] = reference sample file (if DO)
#args[5] = reference legend file (if DO)
#args[6] = covar file

# Input files
# Path to original bamlist
mouse_bamlist <- args[1]

# Path to qtl2 covar file
meta_file <- args[6]
print(meta_file)

# Read covar file
covar <- read.csv(meta_file)

# Use median generation from covar file as the nGen parameter in QUILT

# kill if no generation is included in covar file
stopifnot("gen" %in% colnames(covar))
# kill if there generation covariate isn't numeric or integer
stopifnot(is.numeric(covar$gen) || is.integer(covar$gen))
covar_nGen <- median(covar$gen)

# Chromosome to be analyzed
mouse_chr <- paste0(args[2])

QUILT::QUILT(chr = mouse_chr,
#             regionStart = 50000000, # remove when not testing
#             regionEnd = 55000000, # remove when not testing
#             buffer = 10000, # remove when not testing
             bamlist = mouse_bamlist,
             outputdir = paste0(getwd(), "/"),
             reference_haplotype_file = args[3],
             reference_sample_file = args[4],
             reference_legend_file = args[5],
             nGen = covar_nGen,
             addOptimalHapsToVCF=TRUE,
             record_interim_dosages=TRUE,
             save_prepared_reference=TRUE)
