########################################################################
# Verify that the filtered Sanger VCF really contains only the 7 DO
# founders and biallelic, polymorphic SNPs.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-25
########################################################################

args = commandArgs(trailingOnly = TRUE)

library(VariantAnnotation)

##### VARIABLES #####

base_dir = '/fastscratch/dgatti'

#sanger_file = file.path(base_dir, 'sanger_chr1_do_snps.vcf.gz')
sanger_file = args[1]

##### MAIN #####

# Read in the entire VCF.
vcf = readVcf(sanger_file, 'grcm39')

# How many SNPs do we have?
dim(vcf)

# Verify that the samples are the DO founders. 
colnames(vcf)

# Verify that the variants are all SNPs.
table(info(vcf)$INDEL)

# Verify that the variants are all biallelic & polymorphic in the DO
# founders.
gt = geno(vcf)$GT

# What are the unique calls?
unique(as.vector(gt))

# Are there any duplicated SNP positions?
sum(duplicated(start(rowRanges(vcf))))

# Manually remove rows with missing or heterozygous calls. 
# bcftools didn't do this, but I'm not sure why.
unique_alleles = apply(gt, 1, unique)
remove = grep('0/1', unique_alleles)
vcf = vcf[-remove,]

gt = geno(vcf)$GT
unique_alleles = apply(gt, 1, unique)
remove = grep('\\./\\.', unique_alleles)
vcf = vcf[-remove,]

# What are the unique calls now?
gt = geno(vcf)$GT
unique(as.vector(gt))

# Verify that we only have polymorphic SNPs.
gt = geno(vcf)$GT
num_alt = rowSums(gt == '1/1')
stopifnot(num_alt > 0)

# How many SNPs do we have?
dim(vcf)

# Write out a new VCF with the 7 founders.
# NOTE: writeVcf() uses bgzip, so the suffix is now "bgz".
writeVcf(vcf, filename = sub('\\.gz$', '', sanger_file), index = TRUE)

# Create a data.frame for C57BL/6J and write it out in 23andme format.
# These column names seem to be required by bcftools.
rr = rowRanges(vcf)
df  = data.frame(ID    = names(rr),
                 CHROM = seqnames(rr),
                 POS   = start(rr),
                 AA    = paste0(rr$REF, rr$REF))

write.table(df, file = file.path(base_dir, 'C57BL_6J.tab'), 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

