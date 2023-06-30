################################################################################
# Gather genotypes, markers, & covariates and format for qtl2::read_cross().
# GRCm39.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-28
################################################################################

# library(readxl)
library(VariantAnnotation)
library(qtl2convert)
library(qtl2)

##### VARIABLES #####

# Arguments:
# founder_file: path to founder VCF. 
# sample_file:  path to sample VCF produced by QUILT.
# meta_file:    path to sample metadata file.
# marker_file:  path to QUILT marker file with bp and cM values.
# qtl2_dir:     path to qtl2 output directory where files will be written.
# cross_type:   cross type according to wrt qtl2 preferences

args = commandArgs(trailingOnly = TRUE)

# testing directory
# sample_file_dir_4WC <- "/fastscratch/STITCH_outputDir/work/0e/3600a9168e8fd7e9cb7e8b196c0014"
# sample_file_dir_DO <- "/fastscratch/STITCH_outputDir/work/b4/f2f96658d2d02eeba902552582edd6"

# Founder genotypes and marker positions.
# founder_file_4WC = '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders//chr12_phased_snps.vcf.gz'
# founder_file_DO = '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders/chr12_phased_snps.vcf.gz'
founder_file = args[1]

# Sample genotypes from QUILT.
# sample_file_4WC  = file.path(sample_file_dir_4WC,'pruned.quilt.12.vcf.gz')
# sample_file_DO  = file.path(sample_file_dir_DO,'pruned.quilt.12.vcf.gz')
sample_file = args[2]

# Sample metadata file.
# meta_file_4WC = '/projects/compsci/vmp/lcgbs_ssif/data/GigaMUGA/4WC_covar.csv'
# meta_file_DO = '/projects/compsci/vmp/USERS/widmas/stitch-nf/data/DO_covar.csv'
meta_file = args[3]

# Cross type
# cross_type_4WC = "genail4"
# cross_type_DO = "do"
cross_type = args[4]


# Marker map.
# marker_file_4WC = '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders//chr12_gen_map.txt'
# marker_file_DO = '/projects/compsci/vmp/lcgbs_ssif/data/DO_founders/chr12_gen_map.txt'
marker_file = args[5]

# chromosome
chr = args[6]
# chr = "12"

print(args)


##### MAIN #####

# dir.create(qtl2_dir, showWarnings = FALSE)
print("Read in metadata")
# meta_4WC <- read.csv(meta_file_4WC)
# meta_DO <- read.csv(meta_file_DO)
meta <- read.csv(meta_file)

print("Read in founder genotypes")
# Read in founder genotypes.
founder_vcf = readVcf(founder_file, 'grcm39')
founder_vcf = genotypeCodesToNucleotides(founder_vcf)
founder_gt  = geno(founder_vcf)$GT
founder_gt  = sub('\\|', '', founder_gt)
rownames(founder_gt) = gsub("[^A-Za-z0-9_]", "_", rownames(founder_gt))
head(rownames(founder_gt))

print("Getting founder marker positions")
# Get the marker positions for the founders.
founder_rr = rowRanges(founder_vcf)
names(founder_rr) = gsub("[^A-Za-z0-9_]", "_", names(founder_rr))
print("head(founder_rr)")
head(founder_rr)
rm(founder_vcf)

print("Reading in sample genotypes")
# Read in sample genotypes.
sample_vcf = readVcf(sample_file, 'grcm39')
sample_vcf = genotypeCodesToNucleotides(sample_vcf)
sample_gt  = geno(sample_vcf)$GT
sample_gt  = sub('\\|', '', sample_gt)
stopifnot("id" %in% colnames(meta))
covar_sample_inds <- c()
for(i in 1:length(colnames(sample_gt))){
  ind <- lapply(meta$id, function(x) grep(pattern = x,
                                   x = colnames(sample_gt)[i]))
  covar_sample_inds[i] <- which(ind == 1)
}
stopifnot(length(covar_sample_inds) == length(colnames(sample_gt)))
colnames(sample_gt) <- meta$id[covar_sample_inds]

print("Getting sample marker positions")
# Get the marker positions for the samples.
sample_rr = rowRanges(sample_vcf)
names(sample_rr) = gsub("[^A-Za-z0-9_]", "_", names(sample_rr))
rownames(sample_gt) = gsub("[^A-Za-z0-9_]", "_", rownames(sample_gt))
print("head(sample_rr)")
head(sample_rr)
rm(sample_vcf)

print("Retain SNPs in founder file")
# Retain the SNPs in the founder file.
# We may need to revisit this step later, but I'm being conservative.
common_snps = intersect(rownames(founder_gt), rownames(sample_gt))
founder_gt  = founder_gt[common_snps,]
sample_gt   = sample_gt[common_snps,]
founder_rr  = founder_rr[common_snps,]
sample_rr   = sample_rr[common_snps,]  

# Verify that SNPs line up between founders and samples.
stopifnot(rownames(founder_gt) == rownames(sample_gt))
stopifnot(names(founder_rr)    == names(sample_rr))

# Make heterozygous alleles consistent.
gt_num = t(sample_gt)
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
keep = which(num_geno == 3)
keep = rownames(sample_gt)[keep]

# How many SNPs do we keep?
print("How many alleles do we have?")
length(keep)


print("Filtering founder and sample SNPs")
# Filter the founder and sample SNPs.
founder_gt = founder_gt[keep,]
sample_gt  = sample_gt[keep,]
founder_rr = founder_rr[keep,]
sample_rr  = sample_rr[keep,]

# Verify that SNPs line up between founders and samples.
stopifnot(rownames(founder_gt) == rownames(sample_gt))
stopifnot(names(founder_rr)    == names(sample_rr))

# Merge founders and samples together.
all_gt = base::merge(founder_gt, sample_gt, by = 'row.names',
                     all.x = TRUE, sort = FALSE)
rownames(all_gt) = all_gt$Row.names
all_gt = as.matrix(all_gt[,-1])


print("Get allele codes for each SNP")
# Get the allele codes for each SNP.
ref = as.character(founder_rr$REF)
alt = as.character(unlist(founder_rr$ALT))
alleles = cbind(ref, alt)
rm(ref, alt)

# Encode the genotypes from qtl2.
all_gt = encode_geno(all_gt, alleles)

# Subset the samples to include informative SNPs.
sample_gt = all_gt[,colnames(sample_gt)]

# Write out the founders.
# NOTE: Don't forget to change the order to the standard order!
founder_gt = all_gt[,!colnames(all_gt) %in% colnames(sample_gt)]
colnames(founder_gt) = LETTERS[1:ncol(founder_gt)]
founder_gt = data.frame(marker = rownames(founder_gt), founder_gt) 
write.csv(founder_gt, file = paste0("chr",chr,"_founder_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out sample genotypes.
# keep an eye on the column names here - if there are parsing errors downstream,
# this might be why
sample_gt = data.frame(marker = rownames(sample_gt), 
                       sample_gt, check.names = F)
write.csv(sample_gt, file = paste0("chr",chr,"_sample_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out physical map.
pmap = data.frame(marker = names(founder_rr),
                  chr    = seqnames(founder_rr),
                  pos    = start(founder_rr) * 1e-6)
write.csv(pmap, file = paste0("chr",chr,"_pmap.csv"),
          quote = FALSE, row.names = FALSE)

# Write out genetic map.
markers = read.delim(marker_file, sep = ' ')
markers = subset(markers, position %in% start(founder_rr))
stopifnot(markers$position == start(founder_rr))
gmap = data.frame(marker = names(founder_rr),
                  chr    = seqnames(founder_rr),
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
  stopifnot(colnames(founder_gt)[-1] %in% LETTERS)
  if("sex" %in% colnames(meta)){
    meta$sex = c('female', 'male')[match(meta$sex, c('F', 'M'))]
  } else {
    print("No sex column included")
  }
  write.csv(meta, file = 'covar.csv',
            quote = FALSE, row.names = FALSE)
} else if(cross_type == "do"){
  stopifnot("id" %in% colnames(meta))
  stopifnot("gen" %in% colnames(meta))
  if("sex" %in% colnames(meta)){
    meta$sex = c('female', 'male')[match(meta$sex, c('F', 'M'))]
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

