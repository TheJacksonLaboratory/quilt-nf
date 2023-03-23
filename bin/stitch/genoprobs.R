#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)
library(qtl2fst)

# test_dir
# test_dir <- "/fastscratch/STITCH_outputDir/work/eb/88cc9c53b47f2a073c1485e9f83bcd/"
# setwd(test_dir)

# take arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c("4", "geno4.csv","allele_codes4.csv","pmap4.csv","gmap4.csv","foundergeno4.csv","/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_covar.csv")

# what chromosome?
chr <- args[1]

# read in metadata
metadata <- readr::read_csv(file = args[7], 
                            col_types = c(SampleID = "c",
                                          Sex = "c",
                                          Generation = "n"))

# match sample names to what is encoded in genotype files
geno_header <- c(as.character(readr::read_csv(file = args[2], 
                                              col_names = F, skip = 3, n_max = 1)))[-1]
new_sampleID <- c()
for(i in 1:length(metadata$SampleID)){
  sample_index <- grep(geno_header, pattern = gsub(metadata$SampleID[i], 
                                                   pattern = "-",
                                                   replacement = "."))
  if(length(sample_index) == 0){
    print(paste("sample from metadata not sequenced:", metadata$SampleID[i]))
    new_sampleID[i] <- NA
  } else {
    new_sampleID[i] <- geno_header[sample_index]
  }
}

# Write updated metadata file
updated_DO_metadata <- metadata %>%
  dplyr::mutate(newSampID = new_sampleID) %>%
  dplyr::select(-SampleID) %>%
  dplyr::select(newSampID, everything()) %>%
  dplyr::rename(SampleID = newSampID) %>%
  dplyr::filter(!is.na(SampleID))
write.csv(updated_DO_metadata, file = paste0("chr",chr,"_crossinfo.csv"), quote = F, row.names = F)

# Write control file
qtl2::write_control_file(output_file = paste0("chr",chr,"_control_file.json"),
                         crosstype="do",
                         founder_geno_file=args[6], 
                         founder_geno_transposed=TRUE,
                         gmap_file=args[5],
                         pmap_file=args[4],
                         geno_file=args[2],
                         geno_transposed = TRUE,
                         geno_codes=list(A=1, H=2, B=3), 
                         sex_covar = "Sex",
                         sex_codes=list(F="Female", M="Male"), 
                         covar_file = paste0("chr",chr,"_crossinfo.csv"),
                         crossinfo_covar="Generation",
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

# # Calculate the percent of missing genotypes per sample
# percent_missing <- qtl2::n_missing(cross, "ind", "prop")*100
# missing_genos_df <- data.frame(names(percent_missing), percent_missing) %>%
#   `colnames<-`(c("sample","percent_missing"))

# Calculate genotype probs
# gprobDir <- "/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_seqwell_NovaSeq/geno_probs"
# dir.create(gprobDir, showWarnings = F)
pr <- qtl2::calc_genoprob(cross = cross, 
                          map = cross$pmap, 
                          error_prob = 0.002, 
                          cores = (parallel::detectCores()/2))
# Calculate allele probs
apr <- qtl2::genoprob_to_alleleprob(probs = pr, 
                                    cores = (parallel::detectCores()/2))

# Save objects
save(pr, file = paste0("chr_",chr,"_36_state_probs.RData"))
save(apr, file = paste0("chr_",chr,"_8_state_probs.RData"))




