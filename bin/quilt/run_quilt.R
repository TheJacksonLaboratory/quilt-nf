#!/usr/bin/env Rscript

# Load required packages
library(QUILT)

# Set arguments
args <- commandArgs(trailingOnly = TRUE)

#args[1] = bamlist (from channel)
#args[2] = chromosome number (from channel)
#args[3] = covar file
#args[4] = cross type
#args[5] = reference haplotype file (if DO)
#args[6] = reference sample file (if DO)
#args[7] = reference legend file (if DO)
#args[8] = bin shuffle radius

# Input files
# Path to original bamlist
mouse_bamlist <- args[1]

# Chromosome to be analyzed
mouse_chr <- paste0(args[2])

# Path to qtl2 covar file
meta_file <- args[3]
print(meta_file)

# Read covar file
covar <- read.csv(meta_file)

# Cross type
cross_type <- args[4]

# Reference files
hap <- args[5]
samp <- args[6]
leg <- args[7]
rad <- args[8]

# Take cross type and determine how QUILT should be executed
if(cross_type == "do"){
  
  # DO covar files require generation estimates, so can use this for QUILT
  
  # Use median generation from covar file as the nGen parameter in QUILT
  # kill if no generation is included in covar file
  stopifnot("gen" %in% colnames(covar))
  
  # kill if there generation covariate isn't numeric or integer
  stopifnot(is.numeric(covar$gen) || is.integer(covar$gen))
  stopifnot(length(covar$gen[!is.na(covar$gen)]) > 0)
  covar_nGen <- median(covar$gen[!is.na(covar$gen)])
  
  QUILT::QUILT(chr = mouse_chr,
               #regionStart = 5000000, # remove when not testing
               #regionEnd = 10000000, # remove when not testing
               #buffer = 10000, # remove when not testing
               bamlist = mouse_bamlist,
               outputdir = paste0(getwd(), "/"),
               reference_haplotype_file = hap,
               reference_sample_file = samp,
               reference_legend_file = leg,
               shuffle_bin_radius = as.numeric(rad),
               nGen = covar_nGen,
               addOptimalHapsToVCF=TRUE,
               record_interim_dosages=TRUE,
               save_prepared_reference=TRUE)
  
} else if(cross_type == "bxd"){
  
  # BXD covar files do *not* require generation estimates
  # using roughly 20 generations to fully inbred for now
  covar_nGen <- 20
  
  QUILT::QUILT(chr = mouse_chr,
               #regionStart = 5000000, # remove when not testing
               #regionEnd = 10000000, # remove when not testing
               #buffer = 10000, # remove when not testing
               bamlist = mouse_bamlist,
               outputdir = paste0(getwd(), "/"),
               reference_haplotype_file = hap,
               reference_sample_file = samp,
               reference_legend_file = leg,
               shuffle_bin_radius = as.numeric(rad),
               nGen = covar_nGen,
               addOptimalHapsToVCF=TRUE,
               record_interim_dosages=TRUE,
               save_prepared_reference=TRUE)
  
} else {
  
  print("Cross type specified has no QUILT implementation at this time.")
  
}


