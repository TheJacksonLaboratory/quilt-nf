#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)
library(qtl2fst)


# take arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c("/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_covar.csv", # covars
#           "/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_seqwell_NovaSeq_full/") #sampleDir

# ostem
ostem <- paste0(args[2],"qtl2files")

# read in metadata
metadata <- readr::read_csv(file = args[1], 
                            col_types = c(SampleID = "c",
                                          Sex = "c",
                                          Generation = "n"))

# match sample names to what is encoded in genotype files

# metadata$SampleID <- sample
for_sample_names <- readr::read_csv(file = paste0(ostem,"/geno19.csv"), skip = 3, n_max = 1)
samples <- colnames(for_sample_names)[-1]
matching_samples <- lapply(metadata$SampleID, function(x) samples[grep(x = samples, pattern = x)])
for(i in 1:nrow(metadata)){
  metadata$SampleID[i] <- matching_samples[[i]][1]
}
metadata <- metadata[!is.na(metadata$SampleID),]

# Write updated metadata file
write.csv(metadata, file = paste0(ostem, "/crossinfo.csv"), quote = F, row.names = F)
print("Writing metadata file...")

# Write control file
chr <- seq(1:19)
qtl2::write_control_file(output_file = paste0(ostem,"/lcgbs_control_file.json"),
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
                         covar_file = "crossinfo.csv",
                         crossinfo_covar="Generation",
                         overwrite = T)
print("Writing control file...")

# Load in the cross object
cross <- qtl2::read_cross2(paste0(ostem,"/lcgbs_control_file.json"))
print("Read cross...")

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

# Calculate genotype probs
# pr <- qtl2::calc_genoprob(cross = cross,
#                           map = cross$gmap,
#                           error_prob = 0.002,
#                           quiet = F)
print("Calculating genoprobs...")
fpr <- qtl2fst::calc_genoprob_fst(cross = cross,
                         fbase = "pr", 
                         fdir = paste0(ostem,"/genoprobs"), 
                         error_prob = 0.1,
                         overwrite = T,
                         quiet = F, 
                         cores = (parallel::detectCores()/2))
print("Saving genoprobs..")
save(fpr, cross, file = paste0(ostem,"/lcgbs_36_state_probs.RData"))

# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = fpr, 
                                    cores = (parallel::detectCores()/2))
# Save objects
save(apr, cross, file = paste0(ostem,"/lgcbs_8_state_probs.RData"))


