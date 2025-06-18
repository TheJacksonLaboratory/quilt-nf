#!/usr/bin/env Rscript

################################################################################
# Infer sex of sample from relative coverage on chromosome X.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20250317
################################################################################
library(purrr)
library(dplyr)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
metadata <- args[1]
# metadata <- "/projects/compsci/vmp/USERS/widmas/quilt-nf/data/F1_mut_covar_genail.csv"
covar <- read.csv(metadata, tryLogical = F)

# sample coverage files
gw_cov_files <- list.files(pattern = "GW_coverage")

# chr X coverage files
x_cov_files <- list.files(pattern = "chrX")

# calculate relative coverage of all sites on X vs. all sites in the genome
relativeCoverage <- function(gw, x){
  stopifnot(strsplit(gw,"_")[[1]][[1]] == strsplit(x,"_")[[1]][[1]])
  s_gw <- as.numeric(read.delim(gw,header = F))
  s_x  <- as.numeric(read.delim(x,header = F))
  s <- strsplit(gw,"_")[[1]][[1]]
  rel_cov <- data.frame(s, s_gw, s_x)
  return(rel_cov)
}
relative_coverage_list <- purrr::map2(.x = gw_cov_files,
                                      .y =  x_cov_files,
                                      .f = relativeCoverage)
relative_coverage_df <- Reduce(rbind, relative_coverage_list) %>%
  dplyr::mutate(rel_cov = s_x/s_gw)

# denote sex information
relative_coverage_df <- relative_coverage_df %>%
  dplyr::mutate(inferred_sex = dplyr::if_else(rel_cov > 0.75, "F", "M"))

if(!any(relative_coverage_df$s %in% covar$id)){
  relative_coverage_df$s <- covar$id[unlist(lapply(relative_coverage_df$s, function(x) grep(pattern = x, covar$id)))]
}

# join new sex info to covar
new_covar <- relative_coverage_df %>%
  dplyr::rename(id = s) %>%
  dplyr::right_join(covar,.)

# select variables
sex_checked <- new_covar[,c(colnames(covar),"inferred_sex","rel_cov")] %>%
  dplyr::rename(original_sex = sex,
                sex = inferred_sex)
slim_new_covar <- sex_checked[,c(colnames(covar),"original_sex","rel_cov")]
  
# write the new covar file
write.csv(slim_new_covar, file = "sex_check_covar.csv", row.names = F, quote = F)


