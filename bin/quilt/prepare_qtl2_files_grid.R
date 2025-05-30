################################################################################
# Gather genotypes, markers, & covariates and format for qtl2::read_cross().
# GRCm39.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 2023-09-27
################################################################################

# library(readxl)
library(VariantAnnotation)
library(qtl2convert)
library(qtl2)
library(stringr)


##### VARIABLES #####

# Arguments:
# founder_file: path to founder VCF. 
# sample_file:  path to sample VCF produced by QUILT.
# meta_file:    path to sample metadata file.
# marker_file:  path to QUILT marker file with bp and cM values.
# qtl2_dir:     path to qtl2 output directory where files will be written.
# cross_type:   cross type according to wrt qtl2 preferences

# args = commandArgs(trailingOnly = TRUE)
# 
# # # Founder genotypes and marker positions
# founder_file = args[1]
# # 
# # # Sample genotypes from QUILT.
# sample_file = args[2]
# # 
# # # Sample metadata file.
# meta_file = args[3]
# # 
# # # Cross type
# cross_type = args[4]
# # 
# # # Marker map.
# marker_file = args[5]
# # 
# # # chromosome
# chr = args[6]

# Grid file
grid_file = '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/quilt_1M_physical_grid.csv'


##### TROUBLESHOOTING FILES #####
# Founder genotypes and marker positions
founder_file = '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders/chr19_phased_snps.vcf.gz'
# founder_file = '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/chr19_phased_snps.vcf.gz'

# Sample genotypes from QUILT.
# sample_file = "/projects/reinholdt-lab/DO_ESC/results/quilt/20240224_DO_ESC_downsample/5/2000/quilt_vcfs/quilt.19.vcf.gz"
sample_file = "/projects/compsci/vmp/lcgbs_ssif/results/quilt/20240819_DO_seqwell_grid/0.005/2000/quilt_vcfs/quilt.19.vcf.gz"
# sample_file = "/projects/compsci/vmp/lcgbs_ssif/results/quilt/20240222_4WC_downsampling_gridprobs/3/2000/quilt_vcfs/quilt.19.vcf.gz"

# Sample metadata file.
# meta_file = '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/DO_ESC_covar.csv'
meta_file = '/projects/compsci/vmp/USERS/widmas/quilt-nf/data/DO_covar.csv'
# meta_file = '/projects/compsci/vmp/lcgbs_ssif/data/GigaMUGA/4WC_covar_quilt.csv'

# Cross type
cross_type = 'do'
# cross_type = 'genail4'

# Marker map.
marker_file = '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders/chr19_gen_map.txt'
# marker_file = '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/chr19_gen_map.txt'

# chromosome
chr = "19"




##### MAIN #####

# dir.create(qtl2_dir, showWarnings = FALSE)
print("Read in metadata")
meta <- read.csv(meta_file)

# can't have duplicate sample ids in the metadata for qtl2
meta <- meta[which(!duplicated(meta$id)),]

# can't have any missing values in the covariate file and/or cross info columns
meta <- meta[complete.cases(meta),]

print("Read in founder genotypes")
# Read in founder genotypes.
founder_vcf = readVcf(founder_file, 'grcm39')
founder_vcf = genotypeCodesToNucleotides(founder_vcf)
founder_gt  = geno(founder_vcf)$GT
founder_gt  = sub('\\|', '', founder_gt)
quilt_variants <- nrow(founder_gt)
rownames(founder_gt) = gsub("[^A-Za-z0-9_]", "_", rownames(founder_gt))
head(rownames(founder_gt))

print("Getting founder marker positions")
# Get the marker positions for the founders.
founder_rr = rowRanges(founder_vcf)
names(founder_rr) = gsub("[^A-Za-z0-9_]", "_", names(founder_rr))
# head(founder_rr)
rm(founder_vcf)

print("Reading in sample genotypes")
# Read in sample genotypes.
sample_vcf = readVcf(sample_file, 'grcm39')
sample_vcf = genotypeCodesToNucleotides(sample_vcf)
sample_vcf_info = info(sample_vcf)
rownames(sample_vcf_info) = gsub("[^A-Za-z0-9_]", "_", rownames(sample_vcf_info))

# isolate sample genotypes
sample_gt  = geno(sample_vcf)$GT
sample_gt  = sub('\\|', '', sample_gt)
stopifnot("id" %in% colnames(meta))
covar_sample_inds <- c()

# link the sample IDs 
# assuming that the covar file has some element of the library name in it
covar_sample_ids <- lapply(colnames(sample_gt), function(x){
  # try a symbol split
  first_split <- stringr::str_split(string = x, pattern = "[:punct:]")[[1]]
  sample_search_1 <- lapply(first_split, function(y){
    putative_id <- meta$id[meta$id %in% y]
    if(length(putative_id) > 0){
      return(putative_id)
    }
  })
  result <- unlist(sample_search_1[which(lapply(sample_search_1, function(x) is.null(x)) == FALSE)])
  
  # try a dot split
  dot_split <- stringr::str_split(string = x, pattern = "[.]")[[1]]
  sample_search_2 <- lapply(dot_split, function(y){
    putative_id <- meta$id[meta$id %in% y]
    if(length(putative_id) > 0){
      return(putative_id)
    }
  })
  result2 <- unlist(sample_search_2[which(lapply(sample_search_2, function(x) is.null(x)) == FALSE)])
  
  # try a simple grep
  sample_search_3 <- which(lapply(meta$id, function(y){grep(pattern = paste0("*",y), x = x)}) == 1)
  result3 <- meta$id[sample_search_3]
  
  if(!is.null(result) && length(result) == 1){
    print(result)
    return(result)
  } else if (!is.null(result2) && length(result2) == 1){
    print(result2)
    return(result2)
  } else if (!is.null(result3) && length(result3) == 1){
    print(result3)
    return(result3)
  } else {
    return(NULL)
  }
})

# get rid of samples from quilt vcf not present in the metadata
# remove the null sample names from the covar file ids
quilt_samples_absent_from_meta <- unlist(lapply(covar_sample_ids, is.null))
if(any(quilt_samples_absent_from_meta)){
  covar_sample_ids <- covar_sample_ids[-which(quilt_samples_absent_from_meta)]
  sample_gt <- sample_gt[,-which(quilt_samples_absent_from_meta)]
}

# assign the sample names from the metadata present in both files to the sample genotypes
colnames(sample_gt) <- unlist(covar_sample_ids)

# filter metadata down to samples in the genotype file
meta <- meta[which(meta$id %in% colnames(sample_gt)),]

# at this point, metadata sample size should match the quilt sample size
stopifnot(ncol(sample_gt) == nrow(meta))

print("Getting sample marker positions")
# Get the marker positions for the samples.
sample_rr = rowRanges(sample_vcf)
names(sample_rr) = gsub("[^A-Za-z0-9_]", "_", names(sample_rr))
rownames(sample_gt) = gsub("[^A-Za-z0-9_]", "_", rownames(sample_gt))
rm(sample_vcf)

print("Retain SNPs in founder file")
# Retain the SNPs in the founder file.
# We may need to revisit this step later, but I'm being conservative.
common_snps = intersect(rownames(founder_gt), rownames(sample_gt))
founder_gt  = founder_gt[common_snps,]
sample_gt   = sample_gt[common_snps,]
founder_rr  = founder_rr[common_snps,]
sample_rr   = sample_rr[common_snps,]
sample_vcf_info = sample_vcf_info[rownames(sample_vcf_info) %in% common_snps,]

# Verify that SNPs line up between founders and samples.
stopifnot(rownames(founder_gt) == rownames(sample_gt))
stopifnot(names(founder_rr)    == names(sample_rr))
stopifnot(rownames(founder_gt) == rownames(sample_vcf_info))
stopifnot(names(founder_rr) == rownames(sample_vcf_info))


# load in data from a failed job with this intermediate file in the Nextflow work directory
save(founder_gt, founder_rr, sample_gt, sample_rr, sample_vcf_info, 
     file = paste0("chrom_",chr,"_quilt2qtl2_intermediate.RData"))
# load("/flashscratch/widmas/QUILT_outputDir/work/2a/02e14676a53c74f27cb603c714fe54/chrom_X_quilt2qtl2_intermediate.RData")


# how many sites deviate from HWE?
sample_vcf_info <- data.frame(apply(sample_vcf_info, 2, function(x) unlist(x)))
if(chr != "X"){
  print("Pct of sites that deviate from HWE:")
  paste0(round((table(sample_vcf_info$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")
  # filtering by HWE
  sample_vcf_info <- sample_vcf_info[which(sample_vcf_info$HWE > 0.05),]
} else {
  print("Chr X: Pct of sites that deviate from HWE:")
  print("Skipping this filter")
  paste0(round((table(sample_vcf_info$HWE < 0.05)[[2]]/quilt_variants*100),2),"%")
}


# filtering by info score
lower_info_score <- 0.95
above_threshold_sites <- length(which(sample_vcf_info$INFO_SCORE > lower_info_score))
paste0(signif(above_threshold_sites/length(sample_vcf_info$INFO_SCORE), 4)*100,"% of sites above 0.95 INFO score threshold (",above_threshold_sites,")")

if(above_threshold_sites < 10000){
  print("Fewer than 10,000 sites with info scores > 0.95, setting new threshold and extracting")
  count <- 0
  new_threshold <- lower_info_score
  while (count < 10000) {
    new_threshold <- new_threshold - 0.01  # Increment the threshold
    count <- length(which(sample_vcf_info$INFO_SCORE > new_threshold))
  }
  sample_vcf_info = sample_vcf_info[which(sample_vcf_info$INFO_SCORE > new_threshold),]
  filtered_quilt_variants <- nrow(sample_vcf_info)
} else {
  # This bin is for runs where there were sites above the threshold, but fewer
  print(paste0(above_threshold_sites," info scores above 0.95; keeping just these"))
  sample_vcf_info = sample_vcf_info[which(sample_vcf_info$INFO_SCORE > lower_info_score),]
  filtered_quilt_variants <- nrow(sample_vcf_info)
  new_threshold <- lower_info_score
}

# apply the changes to sample_gt and founder_gt
founder_gt = founder_gt[rownames(founder_gt) %in% rownames(sample_vcf_info),]
sample_gt = sample_gt[rownames(sample_gt) %in% rownames(sample_vcf_info),]
founder_rr = founder_rr[which(names(founder_rr) %in% rownames(sample_vcf_info)),]
sample_rr = sample_rr[which(names(sample_rr) %in% rownames(sample_vcf_info)),]

print("Transposing to grid positions")
grid <- read.csv(grid_file)
grid <- grid[grid$chr == chr,]
grid_gr <- GRanges(seqnames = grid$chr,
                   ranges = IRanges(start = grid$pos*1e6, end = grid$pos*1e6))
grid_nearest <- IRanges::nearest(grid_gr, founder_rr)
grid_nearest <- unique(grid_nearest)

founder_rr_n <- founder_rr[grid_nearest,]
founder_gt_n = founder_gt[rownames(founder_gt) %in% names(founder_rr_n),]
sample_gt_n = sample_gt[rownames(sample_gt) %in% names(founder_rr_n),]
sample_rr_n = sample_rr[which(names(sample_rr) %in% names(founder_rr_n)),]

filtered_quilt_variants <- nrow(sample_gt_n)
resolution_summary <- data.frame(quilt_variants,filtered_quilt_variants, new_threshold)
resolution_summary$chr <- chr
write.csv(resolution_summary, paste0("chr",chr,"_resolution_summary.csv"))

# Make heterozygous alleles consistent.
gt_num = t(sample_gt_n)
gt_num[gt_num == 'CA'] = 'AC'
gt_num[gt_num == 'GA'] = 'AG'
gt_num[gt_num == 'TA'] = 'AT' 
gt_num[gt_num == 'GC'] = 'CG'
gt_num[gt_num == 'TC'] = 'CT'
gt_num[gt_num == 'TG'] = 'GT'


print("Convert genotypes to numerics")
# Convert genotypes to numbers.
# Note that the SNP names will be messed up by the data.frame.
gt_num = data.frame(gt_num, stringsAsFactors = TRUE)
# How many alleles do we have?
print("How many alleles do we have?")
num_geno = sapply(lapply(gt_num, levels), length)
table(num_geno)

# SNPs with one allele just take up space.
# I'm not sure about SNPs with two genotypes. We should look into these and determine
# what's going on. 
# For now, I'm retaining SNPs with all three genotypes.
# Note that I'm using the rownames(sample_gt) because the names of num_geno
# were messed up by data.frame() above.


# # NOTE: when running only 1 sample, this breaks vvv
# if(resolution_summary$filtered_quilt_variants > 35000){
#   keep = which(num_geno == 3)
#   keep = rownames(sample_gt)[keep]
#   
#   # How many SNPs do we keep?
#   print("How many alleles do we have?")
#   length(keep)
# 
#   print("Filtering founder and sample SNPs")
#   # Filter the founder and sample SNPs.
#   founder_gt = founder_gt[keep,]
#   sample_gt  = sample_gt[keep,]
#   founder_rr = founder_rr[keep,]
#   sample_rr  = sample_rr[keep,]
# } else {
#   print("Too few imputed SNPs to restrict to 3 genotypes; keeping what we have")
# }


# Verify that SNPs line up between founders and samples.
stopifnot(rownames(founder_gt_n) == rownames(sample_gt_n))
stopifnot(names(founder_rr_n)    == names(sample_rr_n))

# Merge founders and samples together.
all_gt = base::merge(founder_gt_n, sample_gt_n, by = 'row.names',
                     all.x = TRUE, sort = FALSE)
rownames(all_gt) = all_gt$Row.names
all_gt = as.matrix(all_gt[,-1])


print("Get allele codes for each SNP")
# Get the allele codes for each SNP.
ref = as.character(founder_rr_n$REF)
alt = as.character(unlist(founder_rr_n$ALT))
alleles = cbind(ref, alt)
rm(ref, alt)

# Encode the genotypes from qtl2.
all_gt = encode_geno(all_gt, alleles)

# Subset the samples to include informative SNPs.
sample_gt_n = all_gt[,colnames(sample_gt_n)]

# Write out the founders.
# NOTE: Don't forget to change the order to the standard order!
founder_gt_n = all_gt[,!colnames(all_gt) %in% colnames(sample_gt_n)]
colnames(founder_gt_n) = LETTERS[1:ncol(founder_gt_n)]
founder_gt_n = data.frame(marker = rownames(founder_gt_n), founder_gt_n) 
write.csv(founder_gt_n, file = paste0("chr",chr,"_founder_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out sample genotypes.
# keep an eye on the column names here - if there are parsing errors downstream,
# this might be why
sample_gt_n = data.frame(marker = rownames(sample_gt_n), 
                       sample_gt_n, check.names = F)
write.csv(sample_gt_n, file = paste0("chr",chr,"_sample_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out physical map.
pmap = data.frame(marker = names(founder_rr_n),
                  chr    = seqnames(founder_rr_n),
                  pos    = start(founder_rr_n) * 1e-6)
write.csv(pmap, file = paste0("chr",chr,"_pmap.csv"),
          quote = FALSE, row.names = FALSE)

# Write out genetic map.
markers = read.delim(marker_file, sep = ' ')
markers = subset(markers, position %in% start(founder_rr_n))
stopifnot(markers$position == start(founder_rr_n))
gmap = data.frame(marker = names(founder_rr_n),
                  chr    = seqnames(founder_rr_n),
                  pos    = markers$Genetic_Map.cM.)
write.csv(gmap, file = paste0("chr",chr,"_gmap.csv"),
          quote = FALSE, row.names = FALSE)


# Write out covariates

# Has a cross type been specified?
stopifnot("cross_type" %in% ls(pattern = "cross_type"))

# Does the covar file match the format of the cross type?
if(cross_type == "genail4" || cross_type == "genail8"){
  stopifnot("id" %in% colnames(meta))
  stopifnot("gen" %in% colnames(meta))
  stopifnot(colnames(founder_gt_n)[-1] %in% LETTERS)
  if("sex" %in% colnames(meta)){
    meta$sex[meta$sex == "female" | meta$sex == "f"] <- "F"
    meta$sex[meta$sex == "male" | meta$sex == "m"] <- "M"
  } else {
    print("No sex column included")
  }
  write.csv(meta, file = 'covar.csv',
            quote = FALSE, row.names = FALSE)
} else if(cross_type == "do"){
  stopifnot("id" %in% colnames(meta))
  stopifnot("gen" %in% colnames(meta))
  if("sex" %in% colnames(meta)){
    meta$sex[meta$sex == "female" | meta$sex == "f"] <- "F"
    meta$sex[meta$sex == "male" | meta$sex == "m"] <- "M"
  } else {
    print("No sex column included")
  }
  write.csv(meta, file = 'covar.csv',
            quote = FALSE, row.names = FALSE)
} else {
  "Cross type parsing not included for the specified cross type"
}


# Write out phenotypes.
pheno = data.frame(id  = meta$id,
                   val = rep(1, nrow(meta)))
write.csv(pheno, file = 'pheno.csv',
           quote = FALSE, row.names = FALSE) 

