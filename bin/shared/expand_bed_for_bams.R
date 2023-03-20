#!/usr/bin/env Rscript

# Set arguments
args <- commandArgs(trailingOnly = TRUE)
bed <- args[1]
sample <- args[2]

# read in initial .bed
pile <- read.table(bed, fill = T)

# identify coverage starts
diffs <- c(1,diff(pile[,2]))

# sort the initial .bed from the start sites
start_indices <- which(diffs > 200)
regions <- pile[start_indices,]

# add a 10kb cushion to the starts to create larger window to capture variants
regions$start <- regions$V2 - 10000
regions$end <- regions$V2 + 10000
regions <- regions[,c(1,3,4)]

# write the final .bed file
write.table(regions, file = paste0(sample,"_interval.bed"), row.names = F, col.names = F, quote = F, sep = '\t')