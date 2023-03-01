#!/usr/bin/env Rscript

# Load required packages
library(STITCH)
library(parallel)

# Set arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c("/fastscratch/STITCH_outputDir/work/cc/73ed7fed40b46912004bb84482aa15/STITCH_bamlist.txt",
#           "/fastscratch/STITCH_outputDir/work/cc/73ed7fed40b46912004bb84482aa15/STITCH_15_pos.txt",
#           "4",
#           "15")
# args[1] = bamlist (from channel)
# args[2] = mouse_posfile (from channel)
# args[3] = nFounders
# args[4] = chromosome number (from channel)
# args[5] = generation number if known...hold off on this

# Input files

# Path to original bamlist
mouse_bamlist <- args[1]
temp_bamlist <- read.table(mouse_bamlist)
# mouse_genfile <- paste0(mouse_datadir, "gen.txt")

# Create sample name file in lieu of reheadering .bams
# Avoiding reheadering .bams because initial testing and validation may 
# use simulated or altered data
sample_names <- gsub(unlist(lapply(strsplit(temp_bamlist[,1], 
                                            split = "/"), 
                                   function(x) x[[7]])),
                     pattern = "_dedup.bam", replacement = "")

# Export sample name file
write.table(x = sample_names, 
            file = "sample_names.txt",
            quote = F, 
            row.names = F, 
            col.names = F,
            sep = "\t")
# Define sample names file
sample_names_file <- "sample_names.txt"

# Path to pos file
# pos file is a headerless txt file with four columns:
# 1) chr_ where _ is a number
# 2) position
# 3) ref allele
# 4) alt allele
# mouse_posfile <- "/fastscratch/STITCH_outputDir/work/b7/9f3f7c5a990b6ddc39449bcbb73efd/STITCH_9_pos.txt"

# Read in pos file
temp_posfile <- read.table(args[2])

# Remove duplicate positions
temp_posfile <- temp_posfile[!duplicated(temp_posfile$V2),]

# Remove biallelic multinucleotide variants
SNPs_only <- c(unlist(lapply(strsplit(temp_posfile$V3, ""), function(x) length(x))) == 1)[]
cat(paste0(length(SNPs_only[SNPs_only==F]), 
           " multinucleotide variants removed from position file"))
temp_posfile <- temp_posfile[which(SNPs_only),]

# Write a temporary pos file to be called by STITCH
write.table(x = temp_posfile, 
            file = "temp_pos.txt",
            quote = F, 
            row.names = F, 
            col.names = F,
            sep = "\t")


# Define new pos file
dup_removed_posfile <- "temp_pos.txt"

# Set number of founders
# This should be a param in nextflow
mouse_K <- as.numeric(args[3])

# Number of generations
# This can be a param in nextflow, but maybe not necessary
# mouse_nGen <- args[5]

# Chromosome to be analyzed
mouse_chr <- paste0(args[4])

STITCH::STITCH(tempdir = tempdir(),
               chr = mouse_chr,
               bamlist = mouse_bamlist,
               posfile = dup_removed_posfile,
               outputdir = paste0(getwd(), "/"),
               sampleNames_file = sample_names_file,
               nGen = 20,
               K = mouse_K,
               plotHapSumDuringIterations = FALSE,
               plot_shuffle_haplotype_attempts = FALSE,
               plotAfterImputation = FALSE,
               downsampleToCov = 25)
