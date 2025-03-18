################################################################################
# This script contains code to estimate haplotype blocks in the DO using the
# 36-state genoprobs. It does this by crawling along the probabilities,
# looking for markers which are highly correlated with the current marker.
# The first marker to show low correlation with the current marker is taken
# to be the end of one haplotype block.
# The first function estimates haplotype blocks for ONE SAMPLE on ONE
# CHROMOSOME.
# The second function accepts a qtl2-style genoprobs object and at qtl2-style
# marker map and estimates the haplotype blocks for all samples on all 
# chromosomes.
# The functions are currently SLOW.
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2025-01-16

##########
# Estimate the haplotype blocks for ONE SAMPLE on ONE CHROMOSOME.
# This function crawls along the probability correlation matrix and
# searches for haplotype blocks by finding markers which are highly
# correlated with the first marker in the block. 
# Arguments:
# probs: numeric matrix containing the 36-state genoprobs for
#        one sample on one chromosome. Genotype states in rows
#        and markers in columns. Dimensions must have dimnames.
# mkrs: numeric vector containing the marker positions for the
#       markers in probs. Must be named. Markers must already 
#       be in the same order as probs.
# cor_thr: numeric value indicating the correlation with the proximal 
#          genoprobs at which the end of a hplotype block will be called. 
#          Default = 0.5. 
# alt_prob_thr: numeric value indicating the mean block probability below which
#               alternate genotypes will be listed. Default = 0.9.
#
# Returns: data.frame with one row per haplotype block containing:
#   prox_marker: character string containing the proximal marker.
#   dist_marker: character string containing the distal marker.
#   prox_pos: numeric value containing the proximal marker position.
#   dist_pos: numeric value containing the distal marker position.
#   gt: character string containing the genotype with the highest probability.
#   mean_prob: numeric value containing the mean block probability for the 
#              genotype in "gt".
#   prox_idx: numeric value containing the proximal marker index.
#   dist_idx: numeric value containing the distal marker index.
#   alt_gt: charater string containing other genotype with high probabilities.
#           Each genotype separated by a ":".
#   alt_prob: numeric values containing the mean block probability for the
#             alternate strains listed in "alt_gt". Values separated by ":".
get_blocks_1sample = function(probs, mkrs, cor_thr = 0.5, alt_prob_thr = 0.9) {

  states  = rownames(probs)
  markers = colnames(probs)

  n_blocks = 2000
  blocks = data.frame(prox_marker = rep('', n_blocks),
                      dist_marker = rep('', n_blocks),
                      prox_pos    = rep(0,  n_blocks),
                      dist_pos    = rep(0,  n_blocks),
                      gt          = rep('', n_blocks),
                      mean_prob   = rep(0,  n_blocks),
                      prox_idx    = rep(0,  n_blocks),
                      dist_idx    = rep(0,  n_blocks),
                      alt_gt      = rep(NA, n_blocks),
                      alt_prob    = rep(NA, n_blocks))
  
  block_idx = 1
  prox_idx  = 1
  dist_idx  = 1

  # While the proximal index is less than the number of markers.
  while(prox_idx < length(markers)) {
  
    # Set the search range from the current marker to the end of the
    # chromosome.
    search_range = prox_idx:length(markers)
    
    # Get the correlation of the current marker with all distal markers.
    probs_cor = cor(probs[,prox_idx], probs[,search_range])[1,]
    
    # Find the distal end of the haplotype block by searching for the
    # first marker that has low correlation with the proximal marker.
    wh = which(probs_cor < cor_thr)
    
    # This may happen at the end of the chromosome.
    if(length(wh) == 0) {
    
      dist_idx = length(markers)

    } else {
    
      dist_idx = min(wh) + prox_idx - 2

    } # else

    # If we have a one marker block, skip over it. Otherwise,
    # process the block. 
    if(prox_idx != dist_idx) {
  
      # Set the proximal and distal marker names.
      blocks$prox_marker[block_idx] = markers[prox_idx]
      blocks$dist_marker[block_idx] = markers[dist_idx]
    
      # Get the mean probability for each genotype state.
      rm_probs = rowMeans(probs[,prox_idx:dist_idx])
      
      # Set the called genotype to the state with the highest probability.
      blocks$gt[block_idx] = names(which.max(rm_probs))
      
      # Get the mean probability for the maximal state.
      blocks$mean_prob[block_idx] = max(rm_probs)
    
      # Set the proximal and distal marker indices.
      blocks$prox_idx[block_idx] = prox_idx
      blocks$dist_idx[block_idx] = dist_idx

      # Set the proximan and distal marker positions.
      blocks$prox_pos[block_idx] = mkrs[prox_idx]
      blocks$dist_pos[block_idx] = mkrs[dist_idx]
      
      # Determine whethere we need to include alternate genotypes by
      # seeing whether the mean block probability is over alt_prob_thr.
      rm_probs = sort(rm_probs, decreasing = TRUE)
      # Get the cumulative sum of the genotype probs.
      # We use this to determine whether we should report more than one
      # genotype.
      cs_probs = cumsum(rm_probs)

      num_gt = min(which(cs_probs >= alt_prob_thr))
      if(num_gt > 1) {
      
        # Add alternate genotypes.
        alt_rng = 2:num_gt
        
        blocks$alt_gt[block_idx] = paste(names(rm_probs)[alt_rng], 
                                         collapse = ';')
        blocks$alt_prob[block_idx] = paste(round(rm_probs[alt_rng], digits = 4), 
                                           collapse = ';')
      
      } # if(num_gt > 1)

      # Increment block index.
      block_idx = block_idx + 1
    
    } # if(prox_idx != dist_idx)
    
    # Increment proximal index for the next block.
    prox_idx  = dist_idx  + 1

  } # while(prox_idx < ncol(probs))

  # Trim down the blocks since we started with 1000.
  blocks = subset(blocks, prox_marker != '')
  rownames(blocks) = NULL
  
  return(blocks)

} # get_blocks_1sample()


##########
# This function estimates the haplotype blocks in the DO for all samples
# on all chromosomes.
# Arguments:
# probs: list which is a qtl2-style genoprobs object. Each list element
#        is a 3-dimensional array with samples in rows, 36 genotype states
#        in columns, and markers in slices. All dimensions must have dimnames.
#        The list element names are chromosomes. Marker names must match 
#        those in "map".
# map: list which is a qtl2-style marker map. Each list element is a numeric
#      vector containing marker positions in Mb. The vector names must be
#      marker names. The list element names are chromosomes. Marker names
#      must match those in "probs".
# cor_thr: numeric value indicating the correlation with the proximal 
#          genoprobs at which the end of a hplotype block will be called. 
#          Default = 0.5. 
# alt_prob_thr: numeric value indicating the mean block probability below which
#               alternate genotypes will be listed. Default = 0.9.
# NOTE: I tried to parallelize this with foreach and doParallel, but the
# workers couldn't find the get_blocks_1sample() function even though it
# was in the environment when I called %dopar%.
get_hap_blocks = function(probs, map, cor_thr = 0.5, alt_prob_thr = 0.9) {

  # Verify that probs and map have the same length.
  if(length(probs) != length(map)) {
    stop(paste0('get_hap_blocks: probs (', length(probs), ') and map (',
         length(map), ') must have the same length.'))
  } # if(length(probs) != length(map)

  # Verify that probs and map have the same chromsome names.
  if(!all(names(probs) == names(map))) {
    stop(paste0('get_hap_blocks: probs and map must have the same chromosome ',
         'names in the same order.'))
  } # if(!all(names(probs) == names(map)))

  # Verify that we have a 36-state genoprobs object.
  if(ncol(probs[[1]]) != 36) {
    stop('get_hap_blocks: This function only works with the 36-state genoprobs.')
  } # if(ncol(probs[[1]]) != 36)

  # Verify that the markers are the same between probs and map.
  for(chr in seq_along(probs)) {

    if(!all(names(map[[chr]]) == dimnames(probs[[chr]])[[3]])) {
      stop(paste0('get_hap_blocks: Markers do not match between probs and map',
                  ' on chromosome', names(map)[chr]))
    } # if(!all(names(map[[chr]]) == dimnames(probs[[chr]])[[3]]))

  } # for(chr)

  # Estimate haplotype blocks for all samples on all chromosomes.
  hap_blocks = setNames(vector('list', length(probs)), names(probs))

  for(chr in seq_along(probs)) {

    t1 = proc.time()
    message(paste('CHR:', chr))

    for(sample in 1:nrow(probs[[chr]])) {

      sample_blocks = get_blocks_1sample(probs[[chr]][sample,,], map[[chr]], 
                                         cor_thr, alt_prob_thr)
      sample_blocks = cbind(id  = rownames(probs[[chr]])[sample], 
                            chr = chr, 
                            sample_blocks)

      hap_blocks[[chr]] = rbind(hap_blocks[[chr]], sample_blocks)

    } # for(sample)

    message(proc.time()[3] - t1[3])

  } # for(chr)

  # Combine the haplotype blocks from all chromosomes.
  hap_blocks           = do.call(rbind, hap_blocks)
  rownames(hap_blocks) = NULL
  
  return(hap_blocks)

} # get_hap_blocks()


