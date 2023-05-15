#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)

# test_dir
# test_dir <- "/fastscratch/STITCH_outputDir/work/0e/d3a4eefcb1f765d56edeed34b13d04/"
# setwd(test_dir)

# take arguments
args <- commandArgs(trailingOnly = TRUE)

# what chromosome?
chr <- args[1]
# chr <- "12"

# sample genotypes
sample_genos <- args[2]
# sample_genos <- "chr12_sample_geno.csv"

# founder genotypes
founder_genos <- args[3]
# founder_genos <- "chr12_founder_geno.csv"

# physical map
pmap <- args[4]
# pmap <- "chr12_pmap.csv"

# genetic map
pmap <- args[5]
# gmap <- "chr12_gmap.csv"

# sample metadata
metadata <- args[6]
# metadata <- "covar.csv"


# # read in metadata
# metadata <- readr::read_csv(file = metadata, 
#                             col_types = c(SampleID = "c",
#                                           Sex = "c",
#                                           Generation = "n"))
# 
# # match sample names to what is encoded in genotype files
# geno_header <- c(as.character(readr::read_csv(file = sample_genos, 
#                                               col_names = F, skip = 3, n_max = 1)))[-1]
# new_sampleID <- c()
# for(i in 1:length(metadata$SampleID)){
#   sample_index <- grep(geno_header, pattern = gsub(metadata$SampleID[i], 
#                                                    pattern = "-",
#                                                    replacement = "."))
#   if(length(sample_index) == 0){
#     print(paste("sample from metadata not sequenced:", metadata$SampleID[i]))
#     new_sampleID[i] <- NA
#   } else {
#     new_sampleID[i] <- geno_header[sample_index]
#   }
# }
# 
# # Write updated metadata file
# updated_DO_metadata <- metadata %>%
#   dplyr::mutate(newSampID = new_sampleID) %>%
#   dplyr::select(-SampleID) %>%
#   dplyr::select(newSampID, everything()) %>%
#   dplyr::rename(SampleID = newSampID) %>%
#   dplyr::filter(!is.na(SampleID))
# write.csv(updated_DO_metadata, file = paste0("chr",chr,"_crossinfo.csv"), quote = F, row.names = F)

# Write control file
qtl2::write_control_file(output_file = paste0("chr",chr,"_control_file.json"),
                         crosstype="do",
                         founder_geno_file=founder_genos, 
                         founder_geno_transposed=TRUE,
                         alleles=c(LETTERS[1:8]),
                         gmap_file=gmap,
                         pmap_file=pmap,
                         geno_file=sample_genos,
                         geno_transposed = TRUE,
                         geno_codes=list(A=1, H=2, B=3), 
                         sex_covar = "sex",
                         sex_codes=list(female="female", male="male"), 
                         covar_file = metadata,
                         crossinfo_covar="gen",
                         overwrite = T)

# Load in the cross object
cross <- qtl2::read_cross2(paste0("chr",chr,"_control_file.json"))

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
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = (parallel::detectCores()/2), quiet = F)
# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr, 
                                    cores = (parallel::detectCores()/2))

# Save objects
save(pr, file = paste0("chr_",chr,"_36_state_probs.RData"))
save(apr, file = paste0("chr_",chr,"_8_state_probs.RData"))




