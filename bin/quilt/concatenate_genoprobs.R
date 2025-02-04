#!/usr/bin/env Rscript

################################################################################
# Concatenate genotype probabilities and crosses from QUILT.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250129
################################################################################

library(qtl2)
library(tidyr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)
# test_dir <- "/flashscratch/widmas/QUILT/work/1d/14292a24de96a5e7fd7708ccbbb146/"
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
genoprobs <-  list.files(pattern = "36_state", full.names = T)
genoprobs <- genoprobs[c(1,12:19,2:11,20)]

# cross objects
crosses <- list.files(pattern = "cross", full.names = T)
crosses <- crosses[c(1,12:19,2:11,20)]

# now combining crosses
cross_list <- lapply(crosses, function(x){
  message(x)
  load(x)
  return(cross)
})


# read probs in
probs <- vector(mode = "list", length = length(genoprobs))
names(probs) <- chroms
for(i in 1:length(names(probs))){
  message(genoprobs[i])
  load(genoprobs[i])
  probs[[names(probs)[i]]] <- pr[[1]]
}

# assign attributes
message("Assigning genoprobs attributes...")
attr(probs, "crosstype") <- cross_type
attr(probs, "is_x_chr") <- c(rep(FALSE,19),TRUE)
attr(probs, "alleles") <- cross_list[[1]]$alleles
attr(probs, "alleleprobs") <- FALSE
class(probs) <- c("calc_genoprob", "list")


# combine sample genotypes
message("Combining sample genotypes...")
new_geno <- vector(mode = "list", length = length(cross_list))
names(new_geno) <- chroms
for(i in 1:length(names(new_geno))){
  
  g <- cross_list[[i]]$geno[[1]]
  new_geno[[names(new_geno)[i]]] <- g
}

# combine physical maps
message("Combining physical maps...")
new_pmaps <- vector(mode = "list", length = length(cross_list))
names(new_pmaps) <- chroms
for(i in 1:length(names(new_pmaps))){
  map <- cross_list[[i]]$pmap[[1]]
  new_pmaps[[names(new_pmaps)[i]]] <- map
}

# combine genetic maps
message("Combining genetic maps...")
new_gmaps <- vector(mode = "list", length = length(cross_list))
names(new_gmaps) <- chroms
for(i in 1:length(names(new_gmaps))){
  map <- cross_list[[i]]$gmap[[1]]
  new_gmaps[[names(new_gmaps)[i]]] <- map
}

# combine founder genotypes
message("Combining founder genotypes...")
new_foundergenos <- vector(mode = "list", length = length(cross_list))
names(new_foundergenos) <- chroms
for(i in 1:length(names(new_foundergenos))){
  fg <- cross_list[[i]]$founder_geno[[1]]
  new_foundergenos[[names(new_foundergenos)[i]]] <- fg
}

# skeleton cross to replace with concatenated sample data
cross <- cross_list[[1]]
cross$gmap <- new_gmaps
cross$pmap <- new_pmaps
cross$geno <- new_geno
cross$founder_geno <- new_foundergenos
cross$is_x_chr <- c(rep(FALSE,19),TRUE)

# make allele probs object
message("Generating allele probabilities")
pr <- probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr, quiet = F, cores = parallel::detectCores()/2)

# save everything
message("Saving objects")
saveRDS(object = pr, file = "complete_genoprobs.rds")
saveRDS(object = apr, file = "complete_alleleprobs.rds")
saveRDS(object = cross, file = "complete_cross.rds")



