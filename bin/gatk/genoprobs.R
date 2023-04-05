#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)

# test_dir
# test_dir <- "/fastscratch/STITCH_outputDir/work/ce/abc265dd6791092f92fc8cdfce0e42"
# setwd(test_dir)

# take arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c("plexWell-F04_GT23-02323_ATAGATCC-TTCCTATG_S175_L004", # sample name
#           "/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_covar.csv") # covars

# what sample
sample <- args[1]

# read in metadata
metadata <- readr::read_csv(file = args[2], 
                            col_types = c(SampleID = "c",
                                          Sex = "c",
                                          Generation = "n"))

# match sample names to what is encoded in genotype files
metadata <- metadata[which(lapply(metadata$SampleID, function(x) grep(x = sample, pattern = x))==1),]
metadata$SampleID <- sample
# Write updated metadata file
write.csv(metadata, file = paste0(sample,"_crossinfo.csv"), quote = F, row.names = F)

# Write control file
chr <- seq(1:19)
qtl2::write_control_file(output_file = paste0(sample,"_control_file.json"),
                         crosstype="do",
                         founder_geno_file=paste0("foundergeno",chr,".csv"),
                         founder_geno_transposed=TRUE,
                         gmap_file=paste0("gmap",chr,".csv"),
                         pmap_file=paste0("pmap",chr,".csv"),
                         geno_file=paste0("geno",chr,".csv"),
                         geno_transposed = TRUE,
                         geno_codes=list(A=1, H=2, B=3), 
                         sex_covar = "Sex",
                         sex_codes=list("F"="Female", "M"="Male"), 
                         covar_file = paste0(sample,"_crossinfo.csv"),
                         crossinfo_covar="Generation",
                         overwrite = T)

# Load in the cross object
cross <- qtl2::read_cross2(paste0(sample,"_control_file.json"))

# Drop null markers
cross <- qtl2::drop_nullmarkers(cross)
for(chr in seq_along(cross$founder_geno)) {
  fg <- cross$founder_geno[[chr]]
  g <- cross$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  cross$founder_geno[[chr]] <- fg
  cross$geno[[chr]] <- g
}

# # Calculate the percent of missing genotypes per sample
# percent_missing <- qtl2::n_missing(cross, "ind", "prop")*100
# missing_genos_df <- data.frame(names(percent_missing), percent_missing) %>%
#   `colnames<-`(c("sample","percent_missing"))

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$gmap,
                          error_prob = 0.002, 
                          quiet = F)
# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr, 
                                    cores = (parallel::detectCores()/2))

# Save objects
save(pr, cross, file = paste0(sample,"_36_state_probs.RData"))
save(apr, cross, file = paste0(sample,"_8_state_probs.RData"))


