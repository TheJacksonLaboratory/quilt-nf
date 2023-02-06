#!/usr/bin/env Rscript
library(tictoc)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(qtl2convert)

tictoc::tic()
# Function to encode vcf-style calls in the form of FinalReport files
# Inputs:
# 1) Vector of genotype calls from modified stitch vcf for a single sample
# 2) Sample name corresponding to those genotype calls
# 3) Reference and alternate allele genotypes for each quasimarker (chr+pos)

parseGenos <- function(genos, sample, ref){
  
  # find sample name in genotype and remove it
  binary_geno <- gsub(x = genos, 
                      pattern = paste0(sample,"="), 
                      replacement = "")
  
  # call qtl2-style genotypes from numeric genotype encoding
  binary_geno_df <- data.frame(binary_geno)
  calledGenos <- cbind(binary_geno_df,ref) %>%
    dplyr::mutate(call = dplyr::case_when(binary_geno == "1/1" ~ REF,
                                          binary_geno == "0/0" ~ ALT,
                                          binary_geno %in% c("0/1","1/0") ~ "H",
                                          binary_geno == "./." ~ "-")) %>%
    dplyr::mutate(marker = paste0("st",CHR,"_",POS)) %>%
    dplyr::select(marker, CHR, POS, REF, ALT, call)
  colnames(calledGenos)[6] <- sample
  
  return(calledGenos)
}


# args <- commandArgs(trailingOnly = TRUE)
args <- c("15",
          "/fastscratch/STITCH_outputDir/work/06/57518fb822aac8b90c5c64b3e65180/stitch.15.txt",
          "/fastscratch/STITCH_outputDir/work/08/cb20d2280517dfb1367623e56bfb40/RData/EM.all.15.RData")
# NOTE - need to feed this RData file from somewhere else

# load intermediate data with ancestral genotype calls (inferred founder genotypes)
load(args[3])

# create consensus genotype table
consensusGenos <- cbind(pos, t(eHapsCurrent_tc[, , 1]))
colnames(consensusGenos) <- c(colnames(consensusGenos)[1:4], LETTERS[1:(ncol(consensusGenos)-4)])

# round ancestral genotype priors and calculate which SNPs are segregating among founders
simple_dosage <- consensusGenos %>%
  dplyr::mutate(A = round(A, 1),
                B = round(B, 1),
                C = round(C, 1),
                D = round(D, 1)) %>%
  dplyr::group_by(POS) %>%
  dplyr::summarise(min_dose = min(A, B, C, D),
                   max_dose = max(A, B, C, D),
                   seg = max_dose-min_dose)

# filter the ancestral genotypes to those that segregate among founders
segregating_positions <- simple_dosage %>%
  dplyr::filter(seg > 0.25) %>%
  dplyr::select(POS)

# CREATE ALLELE CODES (qtl2-style)
alleleCodes <- consensusGenos %>%
  dplyr::filter(POS %in% segregating_positions$POS) %>%
  dplyr::select(CHR, POS, REF, ALT) %>%
  dplyr::rename(A = REF,
                B = ALT,
                chr = CHR) %>%
  dplyr::mutate(marker = paste0("st",chr,"_",POS)) %>%
  dplyr::select(marker, chr, `A` ,`B`)


#### SAMPLE GENOTYPES
# read in sample genotypes
# sample_genos <- vroom::vroom(args[2], delim = "\t", n_max = 100000)
sample_genos <- data.table::fread(args[2], 
                                  # nrows = 100000, 
                                  header = TRUE, 
                                  check.names = TRUE)
colnames(sample_genos)[1:4] <- c("CHR","POS","REF","ALT")

# Identify samples
samples <- sample_genos[,5:ncol(sample_genos)]
samples <- as.character(apply(samples, 2, function(x) strsplit(unique(as.matrix(x)), "=")[[1]][1]))

# Reset the column names
colnames(sample_genos) <- c(colnames(sample_genos)[1:4], samples)
sample_genos <- sample_genos[sample_genos$POS %in% segregating_positions$POS]

# remove metadata for parsing
genotypes <- as.matrix(sample_genos[,-c(1:4)])

# Make sample genotypes into a list with each element as a column (sample)
genos_listcols <- lapply(seq_len(ncol(genotypes)), function(i) genotypes[,i])

# Make reference info a list
referenceInfo <- as.data.frame(as.matrix(sample_genos[,c(1:4)]))
referenceInfoList <- lapply(seq_len(length(genos_listcols)), function(x) referenceInfo)

# Supply ref and alt allele calls and encode genotypes in a similar style as FinalReports
parsedGenotypes <- purrr::pmap(.l = list(genos_listcols,
                                    samples,
                                    referenceInfoList),
                               .f = parseGenos)
parsedGenotypes <- suppressMessages(Reduce(dplyr::left_join, parsedGenotypes))
colnames(parsedGenotypes)[1:3] <- c("marker","chr","pos")

# Encode sample genotypes with respect to dynamic founder allele codes
qtl2SampleGenos <- qtl2convert::encode_geno(geno = parsedGenotypes[,-c(1:5)],
                                            allele_codes = alleleCodes[which(alleleCodes$marker %in% parsedGenotypes$marker),c(3:4)])
qtl2SampleGenos <- cbind(parsedGenotypes$marker, data.frame(qtl2SampleGenos))
colnames(qtl2SampleGenos)[1] <- "marker"

# Export sample genotypes in qtl2 format
qtl2convert::write2csv(qtl2SampleGenos, 
                       filename = paste0("geno", args[1], ".csv"), 
                       comment = paste0("stitch-imputed genotypes for chr ", args[1]),
                       overwrite=TRUE)
# Export allele codes in qtl2 format
qtl2convert::write2csv(alleleCodes,
                       filename = paste0("allele_codes",args[1],".csv"), 
                       comment = paste0("Allele codes from STITCH; chromosome", args[1],".csv"),
                       overwrite=TRUE)
tictoc::toc()