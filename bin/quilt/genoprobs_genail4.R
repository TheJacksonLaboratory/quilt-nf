#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)

# test_dir
# setwd(test_dir)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
chrom <- args[1]
# chrom <- "6"

# sample genotypes
sample_genos <- args[2]
# sample_genos <- "chr6_sample_geno.csv"

# founder genotypes
founder_genos <- args[3]
# founder_genos <- "chr6_founder_geno.csv"

# physical map
pmap <- args[4]
# pmap <- "chr6_pmap.csv"

# genetic map
gmap <- args[5]
# gmap <- "chr6_gmap.csv"

# sample metadata
metadata <- args[6]
# metadata <- "covar.csv"
covar <- read.csv(metadata)

# Write control file
# qtl2::write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
#                          crosstype="do",
#                          founder_geno_file=founder_genos,
#                          founder_geno_transposed=TRUE,
#                          alleles=c(LETTERS[1:8]),
#                          gmap_file=gmap,
#                          pmap_file=pmap,
#                          geno_file=sample_genos,
#                          geno_transposed = TRUE,
#                          geno_codes=list(A=1, H=2, B=3),
#                          sex_covar = "sex",
#                          sex_codes=list(female="female", male="male"),
#                          covar_file = metadata,
#                          crossinfo_covar="gen",
#                          overwrite = T)


# Write control file for 4WC mice
write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                   crosstype="genail4",
                   description="4WC_HR",
                   # reading in this cross throws a warning
                   founder_geno_file=founder_genos, 
                   founder_geno_transposed=TRUE,
                   gmap_file=gmap,
                   pmap_file=pmap,
                   geno_file=sample_genos,
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   sex_covar="sex",
                   sex_codes=list(female="female", male="male"),
                   covar_file = metadata,
                   crossinfo_covar = colnames(covar)[!colnames(covar) %in% c("id","sex")],
                   overwrite = T)

# Load in the cross object
cross <- qtl2::read_cross2(paste0("chr",chrom,"_control_file.json"))

# Drop null markers
cross <- qtl2::drop_nullmarkers(cross)
# for(chr in seq_along(cross$founder_geno)) {
#   fg <- cross$founder_geno[[chr]]
#   g <- cross$geno[[chr]]
#   f1 <- colSums(fg==1)/colSums(fg != 0)
#   
#   fg[fg==0] <- NA
#   g[g==0] <- NA
#   
#   fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
#   g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
#   
#   fg[is.na(fg)] <- 0
#   g[is.na(g)] <- 0
#   
#   cross$founder_geno[[chr]] <- fg
#   cross$geno[[chr]] <- g
# }

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = (parallel::detectCores()/2), quiet = F)
# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr, 
                                    cores = (parallel::detectCores()/2))

# Save objects
save(cross, file = paste0("chr_",chrom,"_cross.RData"))
save(pr, file = paste0("chr_",chrom,"_36_state_probs.RData"))
save(apr, file = paste0("chr_",chrom,"_8_state_probs.RData"))







