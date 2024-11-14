library(dplyr)
library(purrr)
library(mmconvert)

# quilt directory
quilt_dir <- "/projects/compsci/vmp/USERS/widmas/quilt-nf"

# chromosome lengths
chr_lengths <- file.path(quilt_dir,"bin/shared/chr_length.txt")
chr_lengths <- read.delim(chr_lengths, header = F)
colnames(chr_lengths) <- c("chr","seq_lengths","isCircular","genome")

# make per-chromosome start and stop
options(scipen = 99999)
thin <- chr_lengths %>%
  dplyr::mutate(start = 3000000) %>%
  dplyr::filter(!chr %in% c("chrM","chrY")) %>%
  dplyr::rename(end = seq_lengths) %>%
  dplyr::select(chr, start, end)
thin$total <- sum(thin$end)

# determine number of markers per chromosome
grid_breakdown_chr <- thin %>%
  dplyr::mutate(pct = end/total,
                n_markers = round(512000*pct,0))
sum(grid_breakdown_chr$n_markers)
grid_breakdown_chr$chr <- gsub(x = grid_breakdown_chr$chr, pattern = "chr", replacement = "")

# get marker grid
grids <- list()
for(i in 1:nrow(thin)){
  start <- grid_breakdown_chr$start[i]
  end <- grid_breakdown_chr$end[i]
  chr <- grid_breakdown_chr$chr[i]
  n_markers <- grid_breakdown_chr$n_markers[i]
  
  pos <- round(seq(from = start, to = end, length.out = n_markers),0)
  print(length(pos))
  grid <- data.frame(chr, pos)
  
  gm_grid <- mmconvert::mmconvert(positions = grid, input_type = "bp")
  
  if(i == 20){
    gm_grid <- gm_grid %>%
      dplyr::select(chr, bp_grcm39, Mbp_grcm39, cM_coxV3_female) %>%
      tidyr::unite(chr:bp_grcm39, col = "marker", remove = F) 
    pmap <- gm_grid %>%
      dplyr::select(chr, marker, Mbp_grcm39) %>%
      dplyr::rename(pos = Mbp_grcm39)
    gmap <- gm_grid %>%
      dplyr::select(chr, marker, cM_coxV3_female) %>%
      dplyr::rename(pos = cM_coxV3_female)
    
    grids[[i]] <- list(pmap, gmap)
  } else {
    gm_grid <- gm_grid %>%
      dplyr::select(chr, bp_grcm39, Mbp_grcm39, cM_coxV3_ave) %>%
      tidyr::unite(chr:bp_grcm39, col = "marker", remove = F) 
    pmap <- gm_grid %>%
      dplyr::select(chr, marker, Mbp_grcm39) %>%
      dplyr::rename(pos = Mbp_grcm39)
    gmap <- gm_grid %>%
      dplyr::select(chr, marker, cM_coxV3_ave) %>%
      dplyr::rename(pos = cM_coxV3_ave)
    
    grids[[i]] <- list(pmap, gmap)
  }
  }
physical_512k_grid <- Reduce(rbind,purrr::transpose(grids)[[1]])
genetic_512k_grid <- Reduce(rbind,purrr::transpose(grids)[[2]])

# write the files
write.csv(physical_512k_grid, file.path(quilt_dir,"data/quilt_512k_physical_grid.csv"), row.names = F, quote = F)
write.csv(genetic_512k_grid, file.path(quilt_dir,"data/quilt_512k_genetic_grid.csv"), row.names = F, quote = F)
