################################################################################
# Use the filtered Sanger VCF, which contains only DO founders and high-quality
# biallelic SNPs, to build a genetic map file for QUILT.
#
# Daniel Gatti & Sam Widmayer
# dan.gatti@jax.org
# samuel.widmayer@jax.org
# 2023-07-11
################################################################################

library(mmconvert)
library(VariantAnnotation)

##### VARIABLES #####

args = commandArgs(trailingOnly = TRUE)

# vcf_dir = 'data/4wc_founders/'
# vcf_file = file.path(vcf_dir, 'sanger_chr1_do_snps.vcf.gz')

vcf_file = args[1]
chr = args[2]

# vcf_file = "data/4wc_founders/chrX_phased_snps.vcf.gz"
# chr = "X"
# outputDir = "data/4wc_founders"

# output_file = file.path(outputDir,paste0("chr",chr,"_gen_map.txt"))

##### MAIN #####

# Read in the filtered Sanger VCF.
vcf = readVcf(vcf_file, 'grcm39')

# Get SNP positions in bp.
pos = data.frame(chr    = seqnames(rowRanges(vcf)),
                 pos    = start(rowRanges(vcf)),
                 marker = names(rowRanges(vcf)))

# Some VCFs have a 'chr' appended to the range; mmconvert doesn't like this
# if(grep(unique(pos$chr), pattern = "chr") == 1){
#   pos$chr <- gsub(x = pos$chr, pattern = "chr", replacement = "")
# } else {
#   print("No 'chr' detected in position df.")
# }

# Use mmconvert to convert from bp to cM.
cm = mmconvert(positions = pos, input_type = 'bp')
if(chr == "X"){
  cm$cM_coxV3_ave = cm$cM_coxV3_female- cm$cM_coxV3_female[1]
} else {
  cm$cM_coxV3_ave = cm$cM_coxV3_ave - cm$cM_coxV3_ave[1]
}


# Create output data frame.
output = data.frame(position = cm$bp_grcm39,
                    COMBINED_rate.cM.Mb. = c(diff(cm$cM_coxV3_ave), 0) / c(diff(cm$Mbp_grcm39), 0) * 1e6,
                    Genetic_Map.cM. = cm$cM_coxV3_ave)
output$COMBINED_rate.cM.Mb.[is.nan(output$COMBINED_rate.cM.Mb.)] = 0

# Write out genetic map file for QUILT.
write.table(output, file = paste0("chr",chr,"_gen_map.txt"), 
            quote = FALSE, row.names = FALSE)
