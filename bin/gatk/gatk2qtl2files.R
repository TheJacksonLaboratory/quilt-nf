#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)
library(qtl2convert)
library(mmconvert)

args <- commandArgs(trailingOnly = TRUE)
# setwd("/fastscratch/STITCH_outputDir/work/97/4bdb6a468db2c6ae6f050e6e800a4b")
# args <- c("18",
#           "/fastscratch/STITCH_outputDir/work/a4/d8d9fa7dcde2d63d00abd1e2b5fdb3/plexWell-F31_GT23-02337_GTTGACAG-TTCCTATG_S189_L004_18_gatk.txt",
#           "/fastscratch/STITCH_outputDir/work/a4/d8d9fa7dcde2d63d00abd1e2b5fdb3/founders_chr18.txt",
#           "8")

#### SAMPLE GENOTYPES
# read in sample genotypes
sample_genos <- data.table::fread(args[2], 
                                  # nrows = 100000, 
                                  header = F, 
                                  check.names = TRUE)
colnames(sample_genos)[1:4] <- c("CHR","POS","REF","ALT")

# read in founder ref genos
reference_genos <- data.table::fread(args[3], 
                                     # nrows = 100000, 
                                     header = F, 
                                     check.names = TRUE)
colnames(reference_genos) <- c("CHR","POS","REF","ALT",LETTERS[1:8])

# attach founder genotypes for calling at the same time
founders_samples <- dplyr::left_join(sample_genos, reference_genos)
null_founders <- which(apply(founders_samples[,(ncol(founders_samples)-8):(ncol(founders_samples))],1,function(x) TRUE %in% is.na(x)))
null_sites <- founders_samples[null_founders,]$POS
founders_samples <- founders_samples[which(!founders_samples$POS %in% null_sites),]

# Identify samples
samples <- founders_samples[,5:ncol(founders_samples)]
samples <- as.character(apply(samples, 2, function(x) strsplit(unique(as.matrix(x)), "=")[[1]][1]))

# Reset the column names
colnames(founders_samples) <- c(colnames(founders_samples)[1:4], samples)
founders_samples <- founders_samples[complete.cases(founders_samples),]

# remove metadata for parsing
genotypes <- as.matrix(founders_samples[,-c(1:4)])

# Make sample genotypes into a list with each element as a column (sample)
genos_listcols <- lapply(seq_len(ncol(genotypes)), function(i) genotypes[,i])

# Make reference info a list
referenceInfo <- as.data.frame(as.matrix(founders_samples[,c(1:4)]))
referenceInfoList <- lapply(seq_len(length(genos_listcols)), function(x) referenceInfo)

# Function to encode vcf-style calls in the form of FinalReport files
# Inputs:
# 1) Vector of genotype calls from modified gatk vcf for a single sample
# 2) Sample name corresponding to those genotype calls
# # 3) Reference and alternate allele genotypes for each quasimarker (chr+pos)
# genos <- genos_listcols[[2]]
# sample <- samples[[2]]
# ref <- referenceInfoList[[2]]
parseGenos <- function(genos, sample, ref){
  
  # find sample name in genotype and remove it
  binary_geno <- gsub(x = genos, 
                      pattern = paste0(sample,"="), 
                      replacement = "")
  # remove phasing :/
  binary_geno <- gsub(x = binary_geno,
                      pattern = "\\|", 
                      replacement = "/")
  
  # call qtl2-style genotypes from numeric genotype encoding
  binary_geno_df <- data.frame(binary_geno)
  calledGenos <- cbind(binary_geno_df,ref) %>%
    dplyr::mutate(call = dplyr::case_when(binary_geno == "0/0" ~ paste0(REF,REF),
                                          binary_geno == "1/1" ~ paste0(ALT,ALT),
                                          binary_geno %in% c("0/1","1/0") ~ paste0(REF,ALT), TRUE ~ "NA")) %>%
    dplyr::mutate(marker = paste0("gatk",CHR,"_",as.numeric(POS))) %>%
    dplyr::select(marker, CHR, POS, REF, ALT, call)
  calledGenos$call[calledGenos$call == "NA"] <- NA
  colnames(calledGenos)[6] <- sample
  
  return(calledGenos)
}

# Supply ref and alt allele calls and encode genotypes in a similar style as FinalReports
options(scipen = 999999999)
parsedGenotypes <- purrr::pmap(.l = list(genos_listcols,
                                         samples,
                                         referenceInfoList),
                               .f = parseGenos)
parsedGenotypes <- suppressMessages(Reduce(dplyr::left_join, parsedGenotypes))
colnames(parsedGenotypes)[1:3] <- c("marker","chr","pos")

# Make allele codes
alleleCodes <- parsedGenotypes %>%
  dplyr::select(marker, chr, REF, ALT) %>%
  dplyr::rename(A = REF,
                B = ALT)

# isolate sample genos
processedSampleGenos <- parsedGenotypes %>%
  dplyr::select(-c(chr, pos, REF, ALT, LETTERS[1:8]))

# isolate founder genos
processedFounderGenos <- parsedGenotypes %>%
  dplyr::select(marker, LETTERS[1:8])


# Encode sample genotypes with respect to dynamic allele codes
qtl2SampleGenos <- qtl2convert::encode_geno(geno = processedSampleGenos[,-1],
                                            allele_codes = alleleCodes[,c(3:4)])
qtl2SampleGenos <- cbind(processedSampleGenos$marker, data.frame(qtl2SampleGenos))
colnames(qtl2SampleGenos) <- colnames(processedSampleGenos)

# Encode founder genotypes with respect to dynamic allele codes
qtl2FounderGenos <- qtl2convert::encode_geno(geno = processedFounderGenos[,-1],
                                            allele_codes = alleleCodes[,c(3:4)])
qtl2FounderGenos <- cbind(processedFounderGenos$marker, data.frame(qtl2FounderGenos))
colnames(qtl2FounderGenos) <- colnames(processedFounderGenos)

# for formatting:
# geno <- rbind(c("C", "G", "C",  "GG", "CG"),
#               c("A", "A", "AT", "TA", "TT"),
#               c("T", "G", NA,   "GT", "TT"))
# codes <- rbind(c("C", "G"), c("A", "T"), c("T", "G"))
# encode_geno(geno, codes)

# Restrict genotypes to sites where founders are not hets
qtl2FounderGenos <- qtl2FounderGenos[which(apply(qtl2FounderGenos, 1, function(x) length(grep("H", x = x))) == 0),]
num_alleles <- apply(qtl2FounderGenos[,-1], 1, function(x) length(summary(as.factor(x))))
seg <- which(num_alleles != 1)
qtl2FounderGenos <- qtl2FounderGenos[seg,]

qtl2SampleGenos <- qtl2SampleGenos[which(qtl2SampleGenos$marker %in% qtl2FounderGenos$marker),]
alleleCodes <- alleleCodes[which(alleleCodes$marker %in% qtl2FounderGenos$marker),]



# Writing qtl2-style files

# Founder Genotypes
qtl2convert::write2csv(qtl2FounderGenos,
                       filename = paste0("foundergeno", args[1], ".csv"), 
                       comment = paste0("gatk-imputed founder genotypes for chr ", args[1]),
                       overwrite=TRUE)

# Sample Genotypes
qtl2convert::write2csv(qtl2SampleGenos, 
                       filename = paste0("geno", args[1], ".csv"), 
                       comment = paste0("gatk-imputed genotypes for chr ", args[1]),
                       overwrite=TRUE)

# Allele Codes
qtl2convert::write2csv(alleleCodes,
                       filename = paste0("allele_codes",args[1],".csv"), 
                       comment = paste0("Allele codes from gatk; chromosome", args[1]),
                       overwrite=TRUE)

# Physical Map
pmap_pos <- as.numeric(unlist(lapply(alleleCodes$marker, function(x) strsplit(x, "_")[[1]][2])))
pmap <- alleleCodes %>%
  dplyr::select(marker, chr) %>%
  dplyr::mutate(pos = pmap_pos/1000000)
qtl2convert::write2csv(pmap,
                       filename = paste0("pmap",args[1],".csv"), 
                       comment = paste0("Physical map from gatk; chromosome",args[1]),
                       overwrite=TRUE)

# Genetic Map
gmap_raw <- mmconvert::mmconvert(positions = pmap %>% dplyr::select(chr, pos, marker), 
                                 input_type = "Mbp")
gmap <- gmap_raw %>%
  dplyr::select(marker, chr, cM_coxV3_ave) %>%
  dplyr::rename(pos = cM_coxV3_ave)
rownames(gmap) <- NULL 
qtl2convert::write2csv(gmap,
                       filename = paste0("gmap",args[1],".csv"), 
                       comment = paste0("Genetic map from gatk; chromosome",args[1]),
                       overwrite=TRUE)
