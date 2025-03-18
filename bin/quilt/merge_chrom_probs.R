#!/usr/bin/env Rscript

################################################################################
# Merge QUILT genotype probabilities from the same chromosome.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250221
################################################################################

library(qtl2)
library(tidyr)
library(dplyr)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
chrom <- args[1]

# genoprob files
pr_files <- list.files(pattern = "pr.rds")

## testing
# test_dir <- "/flashscratch/widmas/QUILT/work/c5/ffbb6a9eff59098539c6f0c930dabe"
# setwd(test_dir)
# chrom <- "1"
# pr_files <- list.files(pattern = "pr.rds")

### RUN ###

# read in probs
probs_list <- vector("list",length(pr_files))
for(i in 1:length(pr_files)){
  message(pr_files[[i]])
  p <- readRDS(pr_files[[i]])
  probs_list[[i]] <- p
}

# Get the depth of the new probs object
array_depth <- sum(unlist(lapply(probs_list, function(x) dim(x)[[3]])))
# Make empty probs object
combined_probs <- array(dim = c(dim(probs_list[[1]])[1], dim(probs_list[[1]])[2], array_depth))
dimnames(combined_probs)[[1]] <- dimnames(probs_list[[1]])[[1]]
dimnames(combined_probs)[[2]] <- dimnames(probs_list[[1]])[[2]]
dimnames(combined_probs)[[3]] <- unlist(lapply(probs_list, function(x) dimnames(x)[[3]]))
# put probs in the right spots with marker names
marker_names <- list()
for(i in 1:length(probs_list)){
  # ind <- new_indices[[i]]
  # message(i)
  # ind <- seq(ind[1], ind[2])
  p <- probs_list[[i]]
  combined_probs[,,dimnames(p)[[3]]] <- p
  # for(j in 1:length(ind)){
  #   combined_probs[,,ind[j]] <- p[,,j]
  # }
  marker_names[[i]] <- dimnames(p)[[3]]
}

# read in marker maps
map_files <- list.files(pattern = "map.rds")
pmap <- unlist(lapply(map_files, readRDS))
pmap <- sort(pmap)

# sort the genoprobs
sorted_combined_probs <- combined_probs[,,names(pmap)]

# write objects
saveRDS(object = sorted_combined_probs, file = paste0("chr",chrom,"_genoprobs.rds"))
saveRDS(object = pmap, file = paste0("chr",chrom,"_pmap.rds"))




