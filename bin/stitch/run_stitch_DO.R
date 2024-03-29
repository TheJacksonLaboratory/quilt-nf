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
# args[5] = reference haplotype file (if DO)
# args[6] = reference sample file (if DO)
# args[7] = reference legend file (if DO)
# args[8] = output file name

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
                                   function(x) x[grep(x, pattern = ".bam")])),
                     pattern = ".bam", replacement = "")

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
            file = paste0("temp_pos_",args[4],".txt"),
            quote = F, 
            row.names = F, 
            col.names = F,
            sep = "\t")


# Define new pos file
dup_removed_posfile <- paste0("temp_pos_",args[4],".txt")

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
               reference_haplotype_file = args[5],
               reference_sample_file = args[6],
               reference_legend_file = args[7],
               nGen = 41,
               K = mouse_K,
               # shuffleHaplotypeIterations = c(5,10,15,20,25,30,40,50,60),
               niterations = 40,
               initial_min_hapProb = 1/mouse_K,
               shuffle_bin_radius = 5000,
               plotHapSumDuringIterations = TRUE,
               plot_shuffle_haplotype_attempts = TRUE,
               plotAfterImputation = TRUE,
               output_haplotype_dosages = TRUE,
               output_filename = args[8])
