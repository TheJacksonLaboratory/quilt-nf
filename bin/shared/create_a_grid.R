library(dplyr)
library(purrr)
library(mmconvert)

# quilt directory
quilt_dir <- "/projects/compsci/vmp/USERS/widmas/quilt-nf"

# cross type
cross_type <- "do"

# grid size
grid_size_M <- 0.25

# chromsome coordinates
ref_maps <- list.files(file.path(quilt_dir,"reference_data",cross_type), pattern = "gen_map", full.names = T)
map_bed <- Reduce(rbind,lapply(ref_maps, function(x){
  start = as.numeric(strsplit(system(command = paste("head -n3", x), intern = T)[[2]]," ")[[1]][[1]])
  stop = as.numeric(strsplit(system(command = paste("tail -n1", x), intern = T)," ")[[1]][[1]])
  chr = gsub(x = strsplit(strsplit(x,"/")[[1]][[length(strsplit(x,"/")[[1]])]],"_")[[1]][[1]], 
             pattern = "chr", 
             replacement = "")
  return(data.frame(chr, start, stop))
}))
map_bed$chr <- factor(map_bed$chr, levels = as.character(c(1:19,"X")))

# get map widths
map_bed <- dplyr::arrange(map_bed, chr) %>%
           dplyr::mutate(width = stop - start)
map_bed$total <- sum(map_bed$width)
map_bed$prop <- map_bed$width/map_bed$total
map_bed$n_markers <- round((grid_size_M*1e6)*map_bed$prop,0)
# if(sum(map_bed$n_markers) < 1e6){
#   rem <- 1e6-sum(map_bed$n_markers)
#   map_bed$n_markers[20] <- map_bed$n_markers[20]+rem
# }

# get marker grid
grids <- list()
for(i in 1:nrow(map_bed)){
  start <- map_bed$start[i]
  end <- map_bed$stop[i]
  chr <- map_bed$chr[i]
  n_markers <- map_bed$n_markers[i]
  
  pos <- seq(from = start, to = end, length.out = n_markers)
  grid <- data.frame(chr, pos)
  
  gm_grid <- mmconvert::mmconvert(positions = grid, input_type = "bp")
  
  if(i == 20){
    gm_grid <- gm_grid %>%
      dplyr::select(chr, bp_grcm39, Mbp_grcm39, cM_coxV3_female) %>%
      dplyr::mutate(bp_grcm39 = as.integer(bp_grcm39)) %>%
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
      dplyr::mutate(bp_grcm39 = as.integer(bp_grcm39)) %>%
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
physical_grid <- Reduce(rbind,purrr::transpose(grids)[[1]])
genetic_grid <- Reduce(rbind,purrr::transpose(grids)[[2]])

# write the files
write.csv(physical_grid, file.path(quilt_dir,"data",paste0("interp_",grid_size_M,"M_physical_grid.csv")), row.names = F, quote = F)
write.csv(genetic_grid, file.path(quilt_dir,"data",paste0("interp_",grid_size_M,"M_genetic_grid.csv")), row.names = F, quote = F)
