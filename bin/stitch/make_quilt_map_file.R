################################################################################
# Use the filtered Sanger VCF, which contains only DO founders and high-quality
# biallelic SNPs, to build a genetic map file for QUILT. Only use Chr 1.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2023-04-23
################################################################################

library(mmconvert)
library(VariantAnnotation)

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

base_dir = '/projects/compsci/vmp/lcgbs_ssif'

vcf_dir = '/fastscratch/widmas'

#vcf_file = file.path(vcf_dir, 'sanger_chr1_do_snps.vcf.gz')
vcf_file = args[1]

output_file = file.path(vcf_dir, 'chr1_gen_map.txt')

##### MAIN #####

# Read in the filtered Sanger VCF.
vcf = readVcf(vcf_file, 'grcm39')

# Get SNP positions in bp.
pos = data.frame(chr    = seqnames(rowRanges(vcf)),
                 pos    = start(rowRanges(vcf)),
                 marker = names(rowRanges(vcf)))

# Use mmconvert to convert from bp to cM.
cm = mmconvert(positions = pos, input_type = 'bp')
cm$cM_coxV3_ave = cm$cM_coxV3_ave - cm$cM_coxV3_ave[1]

# Create output data frame.
output = data.frame(position = cm$bp_grcm39,
                    COMBINED_rate.cM.Mb. = c(diff(cm$cM_coxV3_ave), 0) / c(diff(cm$Mbp_grcm39), 0) * 1e6,
                    Genetic_Map.cM. = cm$cM_coxV3_ave)
output$COMBINED_rate.cM.Mb.[is.nan(output$COMBINED_rate.cM.Mb.)] = 0

# Write out genetic map file for QUILT.
write.table(output, file = output_file, 
            quote = FALSE, row.names = FALSE)
