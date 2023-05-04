################################################################################
# Haplotype Reconstruction on GRCm39 for DO study.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2022-12-06
################################################################################

library(qtl2convert)
library(qtl2)

args = commandArgs(trailingOnly = TRUE)

genome = args[1]

base_dir    = '/projects/compsci/vmp/lcgbs_ssif'
results_dir = file.path(base_dir, 'results', 'quilt')

scratch_dir = '/fastscratch/dgatti'
qtl2_dir    = file.path(scratch_dir, 'qtl2_do_seqwell')

# Haplotype reconstruction.
print('Reading cross')
cross = read_cross2(file.path(qtl2_dir, paste0('seqwell_do_', genome, '.json')))

print(cross)

print('Running HR')
probs = calc_genoprob(cross, map = cross$gmap, cores = 20, quiet = FALSE)
print('Writing 36 state genoprobs')
saveRDS(probs, file = file.path(results_dir, paste0('seqwell_do_genoprobs_', genome, '.rds')))

print('Converting 36 state genoprobs to 8 state allele probs')
aprobs = genoprob_to_alleleprob(probs, cores = 20, quiet = FALSE)
print('Writing allele probs')
saveRDS(aprobs, file = file.path(results_dir, paste0('seqwell_do_alleleprobs_', genome, '.rds')))

