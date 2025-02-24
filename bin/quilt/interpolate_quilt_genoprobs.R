#!/usr/bin/env Rscript

################################################################################
# Interpolate QUILT genotype probabilities to specified gridfile.
#
# Sam Widmayer
# samuel.widmayer@jax.orgs
# 20250221
################################################################################

library(qtl2)
library(tidyr)
library(dplyr)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
chrom <- args[1]

# ranges
start <- args[2]
stop <- args[3]

# qtl2 files
pr_file <- args[4]
cross_file <- args[5]

# grid file
gridfile <- args[6]

# source Gatti interpolation script
source(args[7])

## testing
# test_dir <- "/flashscratch/widmas/QUILT/work/c6/b1acb8cd3ffa61aa16c8329bf4345d"
# setwd(test_dir)
# chrom <- strsplit(list.files(pattern = "cross"),"_")[[1]][[2]]
# start <- "3050143"
# stop <- "9999939"
# # 
# # qtl2 files
# pr_file <- list.files(pattern = "probs")
# cross_file <- list.files(pattern = "cross.RData")
# 
# grid file
# gridfile <- "/projects/compsci/vmp/USERS/widmas/quilt-nf/data/interp_1M_physical_grid.csv"
# source("/projects/compsci/vmp/USERS/widmas/quilt-nf/bin/quilt/interpolate_genoprobs.R")

### RUN ###

# load probs
load(pr_file)

# load cross
load(cross_file)

# load grid
grid <- read.csv(gridfile)
reduced_grid <- grid %>%
  dplyr::filter(chr == chrom,
                pos > min(cross$pmap[[chrom]]),
                pos < max(cross$pmap[[chrom]]))

# make map
grid_map <- reduced_grid$pos
names(grid_map) <- reduced_grid$marker

# interpolate
int_pr <- interpolate_one_chr(pr1 = pr[[chrom]], # original resolution
                              mkr1 = cross$pmap[[chrom]]*1e6, # original map
                              mkr2 = grid_map*1e6) # interpolated probs
saveRDS(object = int_pr, file = paste0("int_",start,"_",stop,"_pr.rds"))
saveRDS(object = grid_map, file = paste0("int_",start,"_",stop,"_map.rds"))
saveRDS(object = cross, file = paste0("original_",start,"_",stop,"_cross.rds"))


