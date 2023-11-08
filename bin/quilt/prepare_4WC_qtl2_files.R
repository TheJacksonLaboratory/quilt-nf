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

args = commandArgs(trailingOnly = TRUE)

# Founder genotypes and marker positions.
# founder_file = '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/chr10_4WC_phased_snps.vcf.gz'
founder_file = args[1]

# Sample genotypes from QUILT.
# sample_file  = '/fastscratch/widmas/pruned.quilt.10.vcf.gz'
sample_file = args[2]

# Sample metadata file.
# meta_file = '/projects/compsci/vmp/lcgbs_ssif/data/GigaMUGA/4WC_covar.csv'
meta_file = args[3]

# Marker map.
# marker_file = '/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders/chr10_gen_map.txt'
marker_file = args[4]

# qtl2 output directory.
# qtl2_dir = getwd()
qtl2_dir = args[5]

# chromosome
chr = args[5]
# chr = 10

print(args)


##### MAIN #####

# dir.create(qtl2_dir, showWarnings = FALSE)

print("Read in founder genotypes")
# Read in founder genotypes.
founder_vcf = readVcf(founder_file, 'grcm39')
founder_vcf = genotypeCodesToNucleotides(founder_vcf)
founder_gt  = geno(founder_vcf)$GT
founder_gt  = sub('\\|', '', founder_gt)
rownames(founder_gt) = sub('\\|', '/', rownames(founder_gt))

print("Getting founder marker positions")
# Get the marker positions for the founders.
founder_rr = rowRanges(founder_vcf)
names(founder_rr) = sub('\\|', '/', names(founder_rr))

rm(founder_vcf)

print("Reading in sample genotypes")
# Read in sample genotypes.
sample_vcf = readVcf(sample_file, 'grcm39')
sample_vcf = genotypeCodesToNucleotides(sample_vcf)
sample_gt  = geno(sample_vcf)$GT
sample_gt  = sub('\\|', '', sample_gt)
rownames(sample_gt)  = gsub(':', '_', rownames(sample_gt))
rownames(sample_gt)  = gsub('/', '_', rownames(sample_gt))

#sample_pos = regexpr('[FM][0-9][0-9]', colnames(sample_gt))
#colnames(sample_gt) = substr(colnames(sample_gt), sample_pos, sample_pos + 2)

print("Getting sample marker positions")
# Get the marker positions for the samples.
sample_rr = rowRanges(sample_vcf)
names(sample_rr) = gsub(':', '_', names(sample_rr))
names(sample_rr) = gsub('/', '_', names(sample_rr))
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
num_geno = sapply(lapply(gt_num, levels), length)
print("How many genotypes have 3 alleles?")
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
length(keep)


print("Filtering founder and sample SNPs to those with three alleles")
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

# Encode the genotypes fro qtl2.
all_gt = encode_geno(all_gt, alleles)

# Subset the samples to include informative SNPs.
sample_gt = all_gt[,-(1:4)]
# TBD: How to filter SNPs with the same SDP.

# Right now, but the chromosome in half....
# keep = 1:floor(nrow(sample_gt) / 2)
# all_gt = all_gt[keep,]
# sample_gt = all_gt[,-(1:8)]

# Read in sample metadata.
meta = read.csv(meta_file)

# Write out the founders.
# NOTE: Don't forget to change the order to the standard order!
founder_gt = all_gt[,1:4]
colnames(founder_gt) = LETTERS[LETTERS %in% colnames(meta)]
founder_gt = data.frame(marker = rownames(founder_gt), founder_gt) 
write.csv(founder_gt, file = paste0("chr",chr,"_founder_geno.csv"),
          quote = FALSE, row.names = FALSE)

# Write out sample genotypes.
sample_gt = data.frame(marker = rownames(sample_gt), sample_gt) 
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

# Match up sample names with those sequenced
# NOTE: we should consider making a sample key
# i.e. sample names according to GigaMUGA; sample names according to lcGBS
# meta = subset(meta, SampleID %in% colnames(sample_gt))
GM_sample_index <- lapply(meta$SampleID, function(x) grep(x, gsub(x = colnames(sample_gt), pattern = "[.]", replacement = "-")))
matched_samples <- lapply(GM_sample_index, function(x){
  if(length(colnames(sample_gt)[x]) == 1){
    colnames(sample_gt)[x]
    } else {
    return(NA)
    }})
meta$SampleID <- unlist(matched_samples)
meta <- meta[which(!is.na(meta$SampleID)),]

# Write out covariates.
# pretty sure this is different for 4WC and CC
# covar = data.frame(id  = meta$SampleID,
#                    sex = meta$Sex,
#                    gen = meta$Generation)
colnames(meta) <- c("id","sex","gen", colnames(meta)[-c(1:3)])
covar <- meta
covar$sex = c('female', 'male')[match(covar$sex, c('F', 'M'))]
write.csv(covar, file = 'covar.csv',
          quote = FALSE, row.names = FALSE)

# Write out phenotypes.
pheno = data.frame(id  = covar$id,
                   val = rep(1, nrow(covar)))
write.csv(pheno, file = 'pheno.csv',
          quote = FALSE, row.names = FALSE)


