#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(purrr)
setwd("/projects/compsci/vmp/USERS/widmas/stitch-nf/")


print("Reading call data")
# gcvf directory
args <- commandArgs(trailingOnly = TRUE)
per_variant_stat_directory <- args[1]
# per_variant_stat_directory <- "data/4WC_ddRADseq_NovaSeq/gvcfs/"

# per-variant stats for each chr
per_variant_stat_files <- list.files(path = per_variant_stat_directory,
                                     pattern = "per_variant_depth_stats.txt")
nomissing_per_variant_stat_files <- list.files(path = per_variant_stat_directory,
                                     pattern = "per_nomissing_variant_depth_stats.txt")

# pull stats
per_variant_coverage <- lapply(per_variant_stat_files, 
                               function(x) data.table::fread(file = paste0(per_variant_stat_directory,x)))
nomissing_per_variant_coverage <- lapply(nomissing_per_variant_stat_files, 
                                         function(x) data.table::fread(file = paste0(per_variant_stat_directory,x)))
biallelic_variants_per_chrom <- lapply(per_variant_coverage, 
                                       function(x) data.frame(unique(x$CHROM), nrow(x))) %>%
  Reduce(rbind,.)
biallelic_variants_per_chrom_nomissing <- lapply(nomissing_per_variant_coverage, 
                                                 function(x) data.frame(unique(x$CHROM), nrow(x))) %>%
  Reduce(rbind,.)

# wrangle stats
colnames(biallelic_variants_per_chrom) <- c("CHROM","n_variants")
biallelic_variants_per_chrom$CHROM <- as.factor(biallelic_variants_per_chrom$CHROM)

colnames(biallelic_variants_per_chrom_nomissing) <- c("CHROM","no_missing_variants")
biallelic_variants_per_chrom_nomissing$CHROM <- as.factor(biallelic_variants_per_chrom_nomissing$CHROM)

# join the two stats
ddRADseq_4WC_variant_counts <- biallelic_variants_per_chrom_nomissing %>%
  dplyr::arrange(CHROM) %>%
  dplyr::full_join(.,biallelic_variants_per_chrom)
genome_wide_variant_counts <- colSums(ddRADseq_4WC_variant_counts[,2:3])

print("Writing number of called variants")
write.csv(ddRADseq_4WC_variant_counts,
          paste0(per_variant_stat_directory,"variant_counts.csv"), 
          quote = F, row.names = F)

# variant counts (millions)
genome_wide_variant_counts/1000000

# calculate average coverage per allele across all samples for no missing sites
print("Calculating per-allele coverage - no missing sites")
sampleVariantAlleleDepth <- function(variant_cov){
  pv_cov <- variant_cov %>%
    tidyr::pivot_longer(-c(CHROM,POS,REF,ALT,DP), names_to = "sample", values_to = "AD")
  print(unique(pv_cov$CHROM))
  pv_cov$sample <- gsub(pv_cov$sample, pattern = ":AD", replacement = "")
  pv_cov$REF_AD <- unlist(lapply(pv_cov$AD, function(x) as.numeric(strsplit(x,",")[[1]][1])))
  pv_cov$ALT_AD <- unlist(lapply(pv_cov$AD, function(x) as.numeric(strsplit(x,",")[[1]][2])))
  pv_cov$totalAD <- pv_cov$REF_AD + pv_cov$ALT_AD
  avg_pv_cov <- pv_cov %>%
    dplyr::group_by(CHROM, POS, REF, ALT, DP) %>%
    dplyr::summarise(mean_REF_AD = mean(REF_AD),
                     mean_ALT_AD = mean(ALT_AD)) %>%
    dplyr::select(CHROM, POS, REF, ALT, DP, mean_REF_AD, mean_ALT_AD)
  return(avg_pv_cov)
}
sample_variant_AD_df <- purrr::map(.x = nomissing_per_variant_coverage,
                                   .f = sampleVariantAlleleDepth) %>%
  Reduce(rbind,.)
sample_variant_AD_df$CHROM <- as.factor(sample_variant_AD_df$CHROM)

per_allele_coverage <- sample_variant_AD_df %>%
  tidyr::pivot_longer(-c(CHROM, POS, REF, ALT, DP), 
                      names_to = "allele_value", 
                      values_to = "AD") %>%
  ggplot(.) + 
  theme_bw() + 
  geom_boxplot(mapping = aes(x = CHROM, y = AD, 
                             colour = allele_value), 
               position = position_dodge(), alpha = 0.3) + 
  ylim(c(0,10)) + 
  ggtitle("Average per-allele coverage - 4WC ddRADseq",
          subtitle = "Only complete sites")
ggsave(per_allele_coverage, filename = paste0(per_variant_stat_directory,"per_allele_coverage.png"), width = 8, height = 8)
  
alt_v_ref <- sample_variant_AD_df %>%
  ggplot(.) + 
  theme_bw() + 
  geom_point(mapping = aes(x = mean_REF_AD, y = mean_ALT_AD), 
             alpha = 0.3) + 
  facet_wrap(.~CHROM, nrow = 4) + 
  xlim(c(0,10)) + 
  ylim(c(0,10)) + 
  ggtitle("Per-allele coverage of ALT vs. REF alleles - 4WC ddRADseq", 
          subtitle = "Only complete sites")
ggsave(alt_v_ref, filename = paste0(per_variant_stat_directory,"alt_v_ref.png"), width = 8, height = 8)


# calculate average coverage per allele across all samples for all sites
print("Calculating per-allele coverage - all called sites")
big_sample_variant_AD_df <- purrr::map(.x = per_variant_coverage,
                                   .f = sampleVariantAlleleDepth) %>%
  Reduce(rbind,.)
big_sample_variant_AD_df$CHROM <- as.factor(big_sample_variant_AD_df$CHROM)

per_allele_coverage <- big_sample_variant_AD_df %>%
  tidyr::pivot_longer(-c(CHROM, POS, REF, ALT, DP), 
                      names_to = "allele_value", 
                      values_to = "AD") %>%
  ggplot(.) + 
  theme_bw() + 
  geom_boxplot(mapping = aes(x = CHROM, y = AD, 
                             colour = allele_value), 
               position = position_dodge(), alpha = 0.3) + 
  ylim(c(0,10)) + 
  ggtitle("Average per-allele coverage - 4WC ddRADseq",
          subtitle = "All calls")
ggsave(per_allele_coverage, filename = paste0(per_variant_stat_directory,"big_per_allele_coverage.png"), width = 8, height = 8)

alt_v_ref <- big_sample_variant_AD_df %>%
  ggplot(.) + 
  theme_bw() + 
  geom_point(mapping = aes(x = mean_REF_AD, y = mean_ALT_AD), 
             alpha = 0.3) + 
  facet_wrap(.~CHROM, nrow = 4) + 
  xlim(c(0,10)) + 
  ylim(c(0,10)) + 
  ggtitle("Per-allele coverage of ALT vs. REF alleles - 4WC ddRADseq", 
          subtitle = "All calls")
ggsave(alt_v_ref, filename = paste0(per_variant_stat_directory,"big_alt_v_ref.png"), width = 8, height = 8)
