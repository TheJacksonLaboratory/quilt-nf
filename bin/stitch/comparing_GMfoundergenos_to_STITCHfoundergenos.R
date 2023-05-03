#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)
library(qtl2)
library(qtl2convert)
library(ggplot2)
library(data.table)

setwd("/fastscratch/widmas/standalone_stitch_testing/")

# SET CHR
chr = 10

# nFounders
nFounders = 8

# FROM STITCH
load(paste0("RData/EM.all.",chr,".10000000.12500000.RData"))

# create consensus genotype table
# each founder column indicates the dosage of the alternate allele
# therefore: -1 = REF/REF, 1 = ALT/ALT
stitchConsensusGenos <- cbind(pos, t(eHapsCurrent_tc[, , 1]))
colnames(stitchConsensusGenos) <- c(colnames(stitchConsensusGenos)[1:4], paste("stitch",LETTERS[1:nFounders], sep = "_"))
rownames(stitchConsensusGenos) <- NULL
stitchConsensusGenos[,-c(1:4)] <- apply(stitchConsensusGenos[,-c(1:4)], 2, function(x) round(x,0))
stitchConsensusGenos <- stitchConsensusGenos[!rowSums(stitchConsensusGenos[,-c(1:4)]) %in% c(0,8),]
stitchConsensusGenos[stitchConsensusGenos==0] <- -1
# pull Sanger founder genotypes
sangerConsensusGenos <- data.frame(data.table::fread(paste0("chr",chr,"_DO_geno_matrix.txt"), 
                                          col.names = c(colnames(stitchConsensusGenos)[1:4], LETTERS[1:nFounders])))
sangerConsensusGenos <- data.frame(sangerConsensusGenos[which(complete.cases(sangerConsensusGenos)),])
sangerConsensusGenos[,-c(1:4)] <- apply(sangerConsensusGenos[,-c(1:4)], 2, function(x) as.numeric(x))
sangerConsensusGenos[,-c(1:4)]$B[1] <- 0

# join founder calls from both sources
all_founder_calls <- dplyr::left_join(stitchConsensusGenos, sangerConsensusGenos)
all_founder_calls <- all_founder_calls[which(complete.cases(all_founder_calls)),]
all_founder_calls$B[1] <- 0

# DENDROGRAM
plot(hclust(dist(cor(sangerConsensusGenos[,-c(1:4,6)]))), main = paste0("chromosome ",chr))

d <- dist(cor(all_founder_calls[,-c(1:4,14)]))
plot(hclust(d), main = paste0("chromosome ",chr))

dm <- as.matrix(d)
founder_model <- Reduce(rbind,apply(dm[-c(1:nFounders),1:nFounders],2, 
      function(x) data.frame(min(x), which(x == min(x))))) %>%
  cbind(.,colnames(dm[-c(1:nFounders),1:nFounders])) %>%
  dplyr::mutate(DO_founder = rownames(.)) %>%
  `colnames<-`(c("d","index","stitch_founder","DO_founder")) %>%
  dplyr::select(-index)
rownames(founder_model) <- NULL
conflicts <- founder_model[grep(founder_model$DO_founder, pattern = "[0-9]"),]
conflict_type <- as.numeric(gsub(x = conflicts$DO_founder, pattern = "[A-Z]", replacement = ""))
if(TRUE %in% conflict_type > 1){
  "same founder tagged three times"
} else {
  missing_stitch_founders <- conflicts$stitch_founder
  missing_DO_founders <- LETTERS[1:nFounders][!LETTERS[1:nFounders] %in% founder_model$DO_founder]
  remaining_founders <- dm[rownames(dm) %in% LETTERS[1:nFounders],colnames(dm) %in% missing_stitch_founders]
  Reduce(rbind,apply(remaining_founders,2,function(x) data.frame(min(x), which(x == min(x))))) %>%
    cbind(.,missing_stitch_founders)
  
  founder_model %>%
    dplyr::filter(!stitch_founder %in% conflicts$stitch_founder) %>%
    dplyr::full_join(.,resolved_founders)
}

# # FROM GIGAMUGA
# GRCm39_founder_genos <- readr::read_csv(paste0("/projects/compsci/vmp/USERS/widmas/DO_Oocyte_Mapping/data/GM_Files/build39/GM_foundergeno",chr,".csv"), skip = 3)
# numeric_GM_founder_genos <- apply(GRCm39_founder_genos[,2:9], c(1,2), function(x) 
#   if(x == "A"){
#     x <- 1
#   } else if(x == "B"){
#     x <- 0
#   } else {
#     x <- NA
#   }
# )
# # DENDROGRAM
# plot(hclust(dist(cor(numeric_GM_founder_genos[complete.cases(numeric_GM_founder_genos),]))), main = "GigaMUGA")

# # attach positions
# grcm39_metadata <- readr::read_csv("/projects/compsci/vmp/USERS/widmas/MUGA_reference_data/data/GigaMUGA/gm_uwisc_v2.csv")
# GRCm39_founder_genos_annot <- GRCm39_founder_genos %>%
#   dplyr::left_join(., grcm39_metadata) %>%
#   dplyr::select(chr, bp_grcm39, everything())
# GRCm39_founder_genos_annot <- GRCm39_founder_genos_annot[,1:11]
# 
# stitch_founder_genos <- example_founder_genos %>%
#   dplyr::mutate(marker = gsub(pattern = "st", replacement = "", marker)) %>%
#   tidyr::separate(marker, into = c("chr", "bp_grcm39"), sep = "_")
# 
# founder_genos_pos <- stitch_founder_genos %>%
#   dplyr::select(chr, bp_grcm39) %>%
#   dplyr::rename(CHROM = chr, 
#                 POS = bp_grcm39)
# 













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
    dplyr::mutate(call = dplyr::case_when(binary_geno == "1/1" ~ paste0(REF,REF),
                                          binary_geno == "0/0" ~ paste0(ALT,ALT),
                                          binary_geno %in% c("0/1","1/0") ~ paste0(REF,ALT), TRUE ~ "NA")) %>%
    dplyr::mutate(marker = paste0("st",CHR,"_",as.numeric(POS))) %>%
    dplyr::select(marker, CHR, POS, REF, ALT, call)
  calledGenos$call[calledGenos$call == "NA"] <- NA
  colnames(calledGenos)[6] <- sample
  
  return(calledGenos)
}

# calls from founder vcf
calls_from_founder_vcf <- readr::read_table(file = "founder_genos_chr6.txt", col_names = FALSE)
colnames(calls_from_founder_vcf) <- c("CHR","POS","REF","ALT",LETTERS[1:8])

# remove metadata for parsing
genotypes <- as.matrix(calls_from_founder_vcf[,-c(1:4)])

# Make sample genotypes into a list with each element as a column (sample)
genos_listcols <- lapply(seq_len(ncol(genotypes)), function(i) genotypes[,i])

# Make reference info a list
referenceInfo <- as.data.frame(as.matrix(calls_from_founder_vcf[,c(1:4)]))
referenceInfoList <- lapply(seq_len(length(genos_listcols)), function(x) referenceInfo)

# Supply ref and alt allele calls and encode genotypes in a similar style as FinalReports
options(scipen = 999999999)
parsedGenotypes <- purrr::pmap(.l = list(genos_listcols,
                                         LETTERS[1:8],
                                         referenceInfoList),
                               .f = parseGenos)
parsedGenotypes <- Reduce(dplyr::left_join, parsedGenotypes)

alleleCodes <- readr::read_csv(paste0("data/DO_seqwell_NovaSeq/qtl2files/allele_codes",chr,".csv"), skip = 3)
encodedFounderGenos <- qtl2convert::encode_geno(geno = parsedGenotypes[,-c(1:5)],
                                            allele_codes = alleleCodes[which(alleleCodes$marker %in% parsedGenotypes$marker),c(3:4)])
encodedFounderGenos <- data.frame(cbind(parsedGenotypes$marker, data.frame(encodedFounderGenos)))
colnames(encodedFounderGenos)[1] <- "marker"

# STITCH founder genos
# strain names unknown
head(example_founder_genos, n = 20)
# from the VCF
# strain names accurate
head(encodedFounderGenos, n = 20)

# Define the bin size
bin_size <- 100
# Create a new column with bin numbers
example_founder_genos$bin <- cut(1:nrow(example_founder_genos), 
                                 breaks = seq(0, nrow(example_founder_genos) + bin_size, bin_size), 
                                 labels = FALSE)
encodedFounderGenos$bin <- cut(1:nrow(encodedFounderGenos), 
                               breaks = seq(0, nrow(encodedFounderGenos) + bin_size, bin_size), 
                               labels = FALSE)

example_founder_genos_binnnested <- example_founder_genos %>%
  # tidyr::pivot_longer(-c(marker,bin), names_to = "founder", values_to = "genotype") %>%
  dplyr::group_by(bin) %>%
  tidyr::nest()

encodedFounderGenos_binnnested <- encodedFounderGenos %>%
  # tidyr::pivot_longer(-c(marker,bin), names_to = "founder", values_to = "genotype") %>%
  dplyr::group_by(bin) %>%
  tidyr::nest()

# Calculate concordant genotypes based on founder VCF and comparing to inferred haplotypes
# stitch_founder_genos <- example_founder_genos_binnnested$data[[996]]
# founder_genos_from_vcf <- encodedFounderGenos_binnnested$data[[996]]
# bin <- encodedFounderGenos_binnnested$bin[[996]]
founderConcordanceCalc <- function(stitch_founder_genos, 
                                   founder_genos_from_vcf,
                                   bin){
  
  # print(bin)
  
  # create founder pair matrix
  founder_pairs <- unique(t(combn(x = c(colnames(stitch_founder_genos)[2:9], colnames(founder_genos_from_vcf)[2:9]), m = 2)))
  
  # comparing stitch ancestral haplotypes to founder vcf
  match <- c()
  for(i in 1:nrow(founder_pairs)){
    # print(i)
    stitch <- stitch_founder_genos[,c("marker",founder_pairs[i,1])]
    vcf <- founder_genos_from_vcf[,c("marker",founder_pairs[i,2])]
    bin_concordance <- dplyr::left_join(stitch, vcf, by = "marker")
    bin_concordance$match <- c(bin_concordance[,2] == bin_concordance[,3])
    perc <- bin_concordance %>%
      dplyr::group_by(match) %>%
      dplyr::count()
    perc$total <- sum(perc$n)
    if(TRUE %in% perc$match){
      match[i] <- perc[which(perc$match == "TRUE"),]$n/perc[which(perc$match == "TRUE"),]$total
    } else {
      match[i] <- 0
    }
    
  }
  out <- data.frame(cbind(founder_pairs, match)) %>%
    dplyr::mutate(bin = bin) %>%
    dplyr::rename(stitch = V1, vcf = V2)
  
  return(out %>%
           dplyr::mutate(match = as.numeric(match)) %>%
           dplyr::group_by(vcf, stitch, bin) %>%
           dplyr::summarise(max = max(match)))
  
}

founder_concordance <- purrr::pmap(.l = list(example_founder_genos_binnnested$data[1:100],
                                             encodedFounderGenos_binnnested$data[1:100],
                                             encodedFounderGenos_binnnested$bin[1:100]), 
                                   .f = founderConcordanceCalc) %>%
  Reduce(dplyr::bind_rows, .)

# save(founder_concordance, file = paste0("founder_concordance_chr",chr,".RData"))

pal <- qtl2::CCcolors
names(pal) <- LETTERS[1:8]
concordance_plot <- founder_concordance %>%
  # dplyr::mutate(diplotype = paste0(stitch, vcf),
  #               bin = as.numeric(bin)) %>%
  ggplot(.) + 
  theme_bw() + 
  geom_smooth(mapping = aes(x = bin, y = mean, colour = vcf)) + 
  scale_colour_manual(values = pal)
concordance_plot
# ggsave(concordance_plot, file = paste0("chr",chr,"_concordance_plot.png"), width = 10, height = 10)

