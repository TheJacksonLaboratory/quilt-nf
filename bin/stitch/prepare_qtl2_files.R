################################################################################
# Gather genotypes, markers, & covariates and format for qtl2::read_cross().
# GRCm39.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-28
################################################################################

library(readxl)
library(VariantAnnotation)
library(qtl2convert)
library(qtl2)

##### VARIABLES #####

# Marker map.
marker_file = '/fastscratch/dgatti/chr1_gen_map.txt'

# Founder genotypes and marker positions.
founder_file = '/fastscratch/dgatti/sanger_chr1_do_snps.vcf.gz'

# Sample genotypes from QUILT.
sample_file  = '/fastscratch/dgatti/quilt_output/quilt.1.vcf.gz'

# Sample metadata file.
meta_file = '/projects/compsci/vmp/lcgbs_ssif/data/DO_tails.xlsx'

# qtl2 output directory.
qtl2_dir = '/fastscratch/dgatti/qtl2'

##### MAIN #####

dir.create(qtl2_dir, showWarnings = FALSE)

# Read in founder genotypes.
founder_vcf = readVcf(founder_file, 'grcm39')
founder_vcf = genotypeCodesToNucleotides(founder_vcf)
founder_gt  = geno(founder_vcf)$GT
founder_gt  = sub('\\|', '', founder_gt)
rownames(founder_gt) = sub('\\|', '/', rownames(founder_gt))

# Get the marker positions for the founders.
founder_rr = rowRanges(founder_vcf)
names(founder_rr) = sub('\\|', '/', names(founder_rr))

rm(founder_vcf)

# Read in sample genotypes.
sample_vcf = readVcf(sample_file, 'grcm39')
sample_vcf = genotypeCodesToNucleotides(sample_vcf)
sample_gt  = geno(sample_vcf)$GT
sample_gt  = sub('\\|', '', sample_gt)
colnames(sample_gt) = gsub('^ddRADseq-|_GT23-[0-9]+_[A-Z]+-[A-Z]+_S[0-9]+$', '', colnames(sample_gt))

# Get the marker positions for the samples.
sample_rr = rowRanges(sample_vcf)

rm(sample_vcf)

# Retain the SNPs in the founder file.
# We may need to revisit this step later, but I'm being conservative.
common_snps = intersect(rownames(founder_gt), rownames(sample_gt))
founder_gt = founder_gt[common_snps,]
sample_gt  = sample_gt[common_snps,]
founder_rr = founder_rr[common_snps,]
sample_rr  = sample_rr[common_snps,]  

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

# Convert genotypes to numbers.
# Note that the SNP names will be messed up by the data.frame.
gt_num = data.frame(gt_num, stringsAsFactors = TRUE)
# How many alleles do we have?
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
length(keep)

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

# Get the allele codes for each SNP.
ref = as.character(founder_rr$REF)
alt = as.character(unlist(founder_rr$ALT))
alleles = cbind(ref, alt)
rm(ref, alt)

# Encode the genotypes fro qtl2.
all_gt = encode_geno(all_gt, alleles)

# Subset the samples to include informative SNPs.
sample_gt = all_gt[,-(1:8)]
colnames(sample_gt) = substr(colnames(sample_gt), 1, 3) 

# Write out the founders.
# NOTE: Don't forget to change the order to the standard order!
founder_gt = all_gt[,1:8]
founder_gt = founder_gt[,c('A_J', 'C57BL_6J', '129S1_SvImJ', 'NOD_ShiLtJ', 
                           'NZO_HlLtJ', 'CAST_EiJ', 'PWK_PhJ', 'WSB_EiJ')]
colnames(founder_gt) = LETTERS[1:8]
founder_gt = data.frame(marker = rownames(founder_gt), founder_gt) 
write.csv(founder_gt, file = file.path(qtl2_dir, 'founder_geno.csv'),
          quote = FALSE, row.names = FALSE)

# Write out sample genotypes.
sample_gt = data.frame(marker = rownames(sample_gt), sample_gt) 
write.csv(sample_gt, file = file.path(qtl2_dir, 'sample_geno.csv'),
          quote = FALSE, row.names = FALSE)

# Write out physical map.
pmap = data.frame(marker = names(founder_rr),
                  chr    = seqnames(founder_rr),
                  pos    = start(founder_rr) * 1e-6)
write.csv(pmap, file = file.path(qtl2_dir, 'pmap.csv'),
          quote = FALSE, row.names = FALSE)

# Write out genetic map.
markers = read.delim(marker_file, sep = ' ')
markers = subset(markers, position %in% start(founder_rr))
stopifnot(markers$position == start(founder_rr))
gmap = data.frame(marker = names(founder_rr),
                  chr    = seqnames(founder_rr),
                  pos    = markers$Genetic_Map.cM.)
write.csv(gmap, file = file.path(qtl2_dir, 'gmap.csv'),
          quote = FALSE, row.names = FALSE)

# Read in sample metadata.
meta = read_xlsx(meta_file)
nc2 = which(nchar(meta$Name) == 2)
meta$Name[nc2] = paste0(substr(meta$Name[nc2], 1, 1), '0', substr(meta$Name[nc2], 2, 2))
meta = subset(meta, Name %in% colnames(sample_gt))

# Write out covariates.
covar = data.frame(id  = meta$Name,
                   sex = substr(meta$Name, 1, 1),
                   gen = meta$Generation)
covar$sex = c('female', 'male')[match(covar$sex, c('F', 'M'))]
write.csv(covar, file = file.path(qtl2_dir, 'covar.csv'),
          quote = FALSE, row.names = FALSE)

# Write out phenotypes.
pheno = data.frame(id  = meta$Name,
                   val = rep(1, nrow(meta)))
write.csv(pheno, file = file.path(qtl2_dir, 'pheno.csv'),
           quote = FALSE, row.names = FALSE) 

# Write out JSON control file.
json = '{
  "description": "ddRADSeq DO",
  "crosstype": "do",
  "sep": ",",
  "na.strings": ["-", "NA"],
  "comment.char": "#",
  "geno": "sample_geno.csv",
  "founder_geno": "founder_geno.csv",
  "gmap": "gmap.csv",
  "pmap": "pmap.csv",
  "pheno": "pheno.csv",
  "covar": "covar.csv",
  "alleles": ["A", "B", "C", "D", "E", "F", "G", "H"],
  "x_chr": "X",
  "genotypes": {
    "A": 1,
    "H": 2,
    "B": 3
  },
  "geno_transposed": true,
  "founder_geno_transposed": true,
  "sex": {
    "covar": "sex",
    "female": "female",
    "male": "male"
  },
  "cross_info": {
    "covar": "gen"
  }
}'

writeLines(text = json, con = file.path(qtl2_dir, 'ddradseq_grcm39.json'))


