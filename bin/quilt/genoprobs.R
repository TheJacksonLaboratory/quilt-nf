#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)

# test_dir
# test_dir <- "/fastscratch/STITCH_outputDir/work/22/7c76bf5b0c3175069a8b8cf3ee3ee1"
# setwd(test_dir)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
chrom <- args[1]
#chrom <- "12"

# sample genotypes
sample_genos <- args[2]
#sample_genos <- "chr12_sample_geno.csv"

# founder genotypes
founder_genos <- args[3]
#founder_genos <- "chr12_founder_geno.csv"

# physical map
pmap <- args[4]
#pmap <- "chr12_pmap.csv"

# genetic map
gmap <- args[5]
#gmap <- "chr12_gmap.csv"

# sample metadata
metadata <- args[6]
#metadata <- "covar.csv"
covar <- read.csv(metadata)

# cross type
cross_type <- args[7]
#cross_type <- "genail4"

# Has a cross type been specified?
stopifnot("cross_type" %in% ls(pattern = "cross_type"))


if(cross_type == "genail4"){
  # write control file for genail4 crosses
  write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                     crosstype="genail4",
                     founder_geno_file=founder_genos,
                     founder_geno_transposed=TRUE,
                     gmap_file=gmap,
                     pmap_file=pmap,
                     geno_file=sample_genos,
                     geno_transposed=TRUE,
                     geno_codes=list(A=1, H=2, B=3),
                     sex_covar="sex",
                     sex_codes=list(female="female", male="male"),
                     covar_file = metadata,
                     crossinfo_covar = colnames(covar)[!colnames(covar) %in% c("id","sex")],
                     overwrite = T)
} else if(cross_type == "do"){
  # Write control file for DO crosses
  qtl2::write_control_file(output_file = paste0("chr",chrom,"_control_file.json"),
                           crosstype="do",
                           founder_geno_file=founder_genos,
                           founder_geno_transposed=TRUE,
                           # alleles=c(LETTERS[1:8]),
                           gmap_file=gmap,
                           pmap_file=pmap,
                           geno_file=sample_genos,
                           geno_transposed = TRUE,
                           geno_codes=list(A=1, H=2, B=3),
                           sex_covar = "sex",
                           sex_codes=list(female="female", male="male"),
                           covar_file = metadata,
                           crossinfo_covar=colnames(covar)[!colnames(covar) %in% "id"],
                           overwrite = T)
} else {
  print("Cross type specified does not have a process to make .json file; ending")
}


# Load in the cross object
cross <- qtl2::read_cross2(paste0("chr",chrom,"_control_file.json"))

# Drop null markers
cross <- qtl2::drop_nullmarkers(cross)

# Calculate genotype probs
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = (parallel::detectCores()/2), quiet = F)
# Estimate genotyping errors
pr_errorlod <- qtl2::calc_errorlod(cross = cross,
                                   probs = pr, 
                                   cores = (parallel::detectCores()/2))

# number of animals with that genotype with an error LOD > 2
genotyping_errors <- rowSums(t(pr_errorlod[[1]]) > 2)
error_genotypes <- names(which(genotyping_errors != 0))
print("How many potential genotyping errors?")
print(paste0(length(error_genotypes), " potential genotyping errors"))

# Remove bad genotypes
corrected_cross <- drop_markers(cross, error_genotypes)
cross <- corrected_cross

# Calculate new genotype probs with reduced map
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap,
                          error_prob = 0.002,
                          cores = (parallel::detectCores()/2), quiet = F)
# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr,
                                    cores = (parallel::detectCores()/2), quiet = F)

# Save objects
save(cross, file = paste0("chr_",chrom,"_cross.RData"))
save(pr, file = paste0("chr_",chrom,"_36_state_probs.RData"))
save(apr, file = paste0("chr_",chrom,"_8_state_probs.RData"))




