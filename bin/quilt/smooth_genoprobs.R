#!/usr/bin/env Rscript

################################################################################
# Count crossovers from QUILT-based haplotype reconstructions using Dan Gatti's
# new correlation-based method.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250210
################################################################################

library(qtl2)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
# test_dir <- "/flashscratch/widmas/QUILT/work/bf/6aaabd153b11c4d95559305af20348"
# setwd(test_dir)

# genotype prob objects
genoprobs <-  args[1]
# genoprobs <- "chr_2_36_state_probs.RData"
load(genoprobs)

# cross objects
cross <- args[2]
# cross <- "chr_2_cross.RData"
load(cross)

# downsample level
downsample_level <- args[3]
# downsample_level <- "0.005"

# bin shuffle param was added at some point but want to keep this flexible
bin_shuffle <- args[4]
# bin_shuffle <- "2000"

# correlation threshold for primary prob call
prob_cor_thresh <- args[5]
# prob_cor_thresh <- "0.5"

# chromosome
chrom <- args[6]
# chrom <- "2"

# source counting function
hap_block_functions <- args[7]
# hap_block_functions <- "/projects/compsci/vmp/USERS/widmas/quilt-nf/bin/quilt/get_hap_blocks.R"
source(hap_block_functions)

# get haplotype blocks one sample at a time
samples = qtl2::ind_ids(cross)
getHapBlocks <- function(ind){

  message(ind)
  
  # probs
  p = pr[[1]][ind,,]
  m = cross$pmap[[1]]
  
  # get haplotype blocks
  haps <- get_blocks_1sample(probs = p, mkrs = m, 
                             cor_thr = prob_cor_thresh)
  haps2 <- haps %>%
    dplyr::mutate(sample = ind)
  return(haps2)
}
hap_block_list <- lapply(samples, getHapBlocks)
hap_blocks_df <- Reduce(rbind,hap_block_list) %>%
  dplyr::mutate(block_size = (as.numeric(dist_pos)-as.numeric(prox_pos)),
                marker_distance = dist_idx-prox_idx)

# identify potentially spurious blocks
flagged_blocks <- hap_blocks_df %>%
  dplyr::filter(marker_distance < 5) %>%
  dplyr::filter(block_size < 0.01)

# grab marker names
idx_list <- vector("list", nrow(flagged_blocks))
for(i in 1:nrow(flagged_blocks)){
  idx_df <- data.frame(seq(flagged_blocks[i,]$prox_idx,flagged_blocks[i,]$dist_idx),
             flagged_blocks[i,]$sample)
  colnames(idx_df) <- c("idx","sample")
  idx_list[[i]] <- idx_df
}
idx_samples_df <- Reduce(rbind,idx_list)
flagged_markers <- sort(unique(idx_samples_df$idx))
flagged_markers <- names(cross$pmap[[1]][flagged_markers])

# remove the spurious markers
cross <- qtl2::drop_markers(cross, markers = flagged_markers)

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = (parallel::detectCores()/2), quiet = F)

# Save objects
save(cross, file = paste0("chr_",chrom,"_cross_smooth.RData"))
save(pr, file = paste0("chr_",chrom,"_36_state_probs_smooth.RData"))

