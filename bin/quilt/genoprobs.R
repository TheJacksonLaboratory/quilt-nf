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

# # temporary for integrating
# load("results/quilt/20231120_DO_seqwell_downsampled/3/geno_probs/chr_19_cross.RData")
# load("results/quilt/20231120_DO_seqwell_downsampled/3/geno_probs/chr_19_36_state_probs.RData")
# load("results/quilt/20231120_DO_seqwell_downsampled/3/geno_probs/chr_19_8_state_probs.RData")

print("Counting Crossovers.")
# count crossovers
m_quilt <- qtl2::maxmarg(probs = pr,
                         cores = c(parallel::detectCores()/2))
crossovers <- qtl2::count_xo(geno = m_quilt, 
                             quiet = F, 
                             cores = c(parallel::detectCores()/2))

# find recombination breakpoints
crossover_locs <- qtl2::locate_xo(geno = m_quilt,
                                  map = cross$pmap, 
                                  cores = c(parallel::detectCores()/2))[[1]]

# measuring recombination blocks
xo_block_sizes <- lapply(1:length(crossover_locs), function(ind){
  # print(ind)
  sample <- names(crossover_locs[ind])
  if(length(crossover_locs[[ind]]) == 0){
    print("No crossovers detected on individual; skipping")
  } else {
    xo_lengths <- list()
    for(xo in 1:length(crossover_locs[[ind]])){
      if(xo == 1){
        df <- data.frame(sample,crossover_locs[[ind]][xo] - cross$pmap[[1]][[1]])
        df$block_start <- cross$pmap[[1]][[1]]
      } else {
        df <- data.frame(sample,crossover_locs[[ind]][xo] - crossover_locs[[ind]][xo-1])
        df$block_start <- crossover_locs[[ind]][xo-1]
      }
      df$block_end <- crossover_locs[[ind]][xo]
      colnames(df)[1:2] <- c("sample","block_size")
      xo_lengths[[xo]] <- df
    }
    return(Reduce(dplyr::bind_rows,xo_lengths))
  }
})

# find recombination blocks < a certain size
find_spurious_xo_markers <- purrr::map2(.x = crossover_locs, .y = xo_block_sizes, .f = function(x,y){
  small_blocks = x[which(y$block_size < 0.001)]
  if(length(small_blocks) == 0){
    print("No crossovers with length below threshold")
  } else {
    markers_to_remove <- qtl2::find_marker(map = cross$pmap, 
                                           chr = chrom,
                                           pos = small_blocks)
  }
  return(markers_to_remove)})
names(find_spurious_xo_markers) <- NULL
spurious_xo_markers <- unlist(find_spurious_xo_markers)

# how many markers contribute this?
print(paste("Unique spurious markers:", length(unique(spurious_xo_markers))))
print(paste("Total spurious markers:", length(spurious_xo_markers)))

# remove these markers from the requisite cross object elements
cross <- drop_markers(cross = cross, unique(spurious_xo_markers))

# Calculate new genotype probs with reduced map AGAIN
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