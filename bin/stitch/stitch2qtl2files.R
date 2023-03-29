#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(purrr)
library(qtl2)
library(qtl2convert)
library(mmconvert)

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
    dplyr::mutate(call = dplyr::case_when(binary_geno == "0/0" ~ paste0(REF,REF),
                                          binary_geno == "1/1" ~ paste0(ALT,ALT),
                                          binary_geno %in% c("0/1","1/0") ~ paste0(REF,ALT), TRUE ~ "NA")) %>%
    dplyr::mutate(marker = paste0("st",CHR,"_",as.numeric(POS))) %>%
    dplyr::select(marker, CHR, POS, REF, ALT, call)
  calledGenos$call[calledGenos$call == "NA"] <- NA
  colnames(calledGenos)[6] <- sample
  
  return(calledGenos)
}

# Function to convert allele dosages among founders from STITCH to
# letter code genotypes
# Inputs:
# 1) Data frame with columns CHR, POS, REF, ALT, and founder name (typically A, B, C, D, etc)
# founder name column has dosage values

founderDoseToAllele <- function(founder_dose){
  parent <- colnames(founder_dose)[5]
  colnames(founder_dose)[5] <- "dose"
  founder_call<- founder_dose %>%
    dplyr::mutate(call = dplyr::case_when(dose < 0.25 ~ as.character(REF),
                                          dose > 0.75 ~ as.character(ALT),
                                          TRUE ~ "H"))
  colnames(founder_call)[6] <- parent
  founder_call <- founder_call %>% dplyr::select(-dose)
  return(founder_call)
  
}

args <- commandArgs(trailingOnly = TRUE)
# args <- c("17",
#           "/fastscratch/STITCH_outputDir/work/62/8c80c56fe94ad5b6c15a749f18a4e3/stitch.17.txt",
#           "data/DO_seqwell_NovaSeq/stitch_vcfs/RData/EM.all.17.RData",
#           "8")

# load intermediate data with ancestral genotype calls (inferred founder genotypes)
load(args[3])

# encode estimated number of founders
nFounders <- as.numeric(args[4])

# create consensus genotype table
# each founder column indicates the dosage of the alternate allele
# therefore: 0 = REF/REF, 1 = ALT/ALT
consensusGenos <- cbind(pos, t(eHapsCurrent_tc[, , 1]))
colnames(consensusGenos) <- c(colnames(consensusGenos)[1:4], LETTERS[1:nFounders])

# round ancestral genotype priors and calculate which SNPs are segregating among founders
roundedConsensusGenos <- apply(consensusGenos[-c(1:4)], 2, function(x) round(x, digits = 1))
seg <- apply(roundedConsensusGenos, 1, function(x) max(x)-min(x))

simple_dosage <- cbind(consensusGenos[,c(1:4)], roundedConsensusGenos) %>%
  dplyr::mutate(seg = seg)

# filter the ancestral genotypes to those that segregate among founders
segregating_positions <- simple_dosage %>%
  dplyr::filter(seg > (1/nFounders)) %>%
  dplyr::select(POS)

# Create allele codes (qtl2-style)
alleleCodes <- consensusGenos %>%
  dplyr::filter(POS %in% segregating_positions$POS) %>%
  dplyr::select(CHR, POS, REF, ALT) %>%
  dplyr::rename(A = REF,
                B = ALT,
                chr = CHR) %>%
  dplyr::mutate(marker = paste0("st",chr,"_",POS)) %>%
  dplyr::select(marker, chr, `A` ,`B`)

# Create founder consensus genos
# genotype is calculated from dose of alt allele (quantitative), so need 
# a way to call alleles

# filter to segregating sites
raw_founder_calls <- consensusGenos %>%
  dplyr::filter(POS %in% segregating_positions$POS)

# make a list of data frames with allele doses for each founder
founder_calls_list <- lapply(seq_len(nFounders), function(x) raw_founder_calls[,c(1:4,x+4)])

# call genotypes from doses
consensusGenosLetters <- suppressMessages(purrr::map(.x = founder_calls_list, 
                                                     .f = founderDoseToAllele) %>%
                                            Reduce(left_join,.)) %>%
  dplyr::mutate(marker = paste0("st",CHR,"_",POS)) %>%
  dplyr::select(-REF, -ALT) %>%
  dplyr::select(marker, CHR, POS, everything()) %>%
  dplyr::rename(chr = CHR,
                pos = POS)

# recode genotypes using the allele codes
recodedConsensusGenos <- qtl2convert::encode_geno(geno = consensusGenosLetters[c(((ncol(consensusGenosLetters)-nFounders)+1):ncol(consensusGenosLetters))], 
                         allele_codes = alleleCodes[which(alleleCodes$marker %in% consensusGenosLetters$marker),c(3:4)])
recodedConsensusGenos <- cbind(raw_founder_calls[,c(1:4)], data.frame(recodedConsensusGenos)) %>%
  dplyr::select(-REF, -ALT) %>%
  dplyr::rename(chr = CHR,
                pos = POS) %>%
  dplyr::mutate(marker = paste0("st",chr,"_",pos)) %>%
  dplyr::select(-chr, -pos) %>%
  dplyr::select(marker, everything())


#### SAMPLE GENOTYPES
# read in sample genotypes
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
options(scipen = 999999999)
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
# for formatting:
# geno <- rbind(c("C", "G", "C",  "GG", "CG"),
#               c("A", "A", "AT", "TA", "TT"),
#               c("T", "G", NA,   "GT", "TT"))
# codes <- rbind(c("C", "G"), c("A", "T"), c("T", "G"))
# encode_geno(geno, codes)


# Writing qtl2-style files

# Founder Genotypes
qtl2convert::write2csv(recodedConsensusGenos,
                       filename = paste0("foundergeno", args[1], ".csv"), 
                       comment = paste0("stitch-imputed founder genotypes for chr ", args[1]),
                       overwrite=TRUE)

# Sample Genotypes
qtl2convert::write2csv(qtl2SampleGenos, 
                       filename = paste0("geno", args[1], ".csv"), 
                       comment = paste0("stitch-imputed genotypes for chr ", args[1]),
                       overwrite=TRUE)

# Allele Codes
qtl2convert::write2csv(alleleCodes,
                       filename = paste0("allele_codes",args[1],".csv"), 
                       comment = paste0("Allele codes from STITCH; chromosome", args[1],".csv"),
                       overwrite=TRUE)

# Physical Map
pmap_pos <- as.numeric(unlist(lapply(alleleCodes$marker, function(x) strsplit(x, "_")[[1]][2])))
pmap <- alleleCodes %>%
  dplyr::select(marker, chr) %>%
  dplyr::mutate(pos = pmap_pos/1000000)
qtl2convert::write2csv(pmap,
                       filename = paste0("pmap",args[1],".csv"), 
                       comment = paste0("Physical map from STITCH; chromosome", args[1],".csv"),
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
                       comment = paste0("Genetic map from STITCH; chromosome", args[1],".csv"),
                       overwrite=TRUE)
