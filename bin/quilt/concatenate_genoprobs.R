#!/usr/bin/env Rscript

################################################################################
# Concatenate genotype probabilities and physical maps from QUILT.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250221
################################################################################

library(qtl2)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
# test_dir <- "/flashscratch/widmas/QUILT/work/71/9f293f2781e9ad38bc9a497f0d2111"
# cross_type = "do"
# setwd(test_dir)

# cross type
cross_type <- args[1]
if(cross_type == "het3" | cross_type == "cc" | cross_type == "genail4"){
  cross_type <- "genail4"
} else if(cross_type == "bxd"){
  cross_type <- "risib"
} else if(cross_type == "do"){
  cross_type <- "do"
} else {
  "No clue of cross type!"
}

# chromosomes
chroms <- c(as.character(seq(1:19)),"X")

# genotype prob objects
genoprobs <-  paste0("chr",chroms,"_genoprobs.rds")

# read probs in
probs <- vector(mode = "list", length = length(genoprobs))
names(probs) <- chroms
for(i in 1:length(names(probs))){
  message(genoprobs[i])
  pr <- readRDS(genoprobs[i])
  probs[[names(probs)[i]]] <- pr
}

# assign attributes
message("Assigning genoprobs attributes...")
attr(probs, "crosstype") <- cross_type
attr(probs, "is_x_chr") <- c(rep(FALSE,19),TRUE)
attr(probs, "alleleprobs") <- FALSE
class(probs) <- c("calc_genoprob", "list")
attr(probs, "alleles") <- unique(unlist(lapply(dimnames(probs[[1]])[[2]], 
                                               function(x) strsplit(x, split = "")[[1]])))

# combine physical maps
message("Combining physical maps...")
pmaps <-  paste0("chr",chroms,"_pmap.rds")
new_pmaps <- vector(mode = "list", length = length(pmaps))
names(new_pmaps) <- chroms
for(i in 1:length(names(new_pmaps))){
  map <- readRDS(pmaps[i])
  new_pmaps[[names(new_pmaps)[i]]] <- map
}

# make allele probs object
message("Generating allele probabilities")
pr <- probs
rm(probs)
apr <- qtl2::genoprob_to_alleleprob(probs = pr, quiet = F, cores = parallel::detectCores()/1.2)

# save everything
message("Saving objects")
saveRDS(object = pr, file = "complete_genoprobs.rds")
saveRDS(object = apr, file = "complete_alleleprobs.rds")
saveRDS(object = new_pmaps, file = "complete_pmap.rds")



