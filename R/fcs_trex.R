fcs_trex <- function(fcs_join_obj, compare_list, reduction = c("UMAP","tSNE"), outdir,
                     point_alpha = 0.05, neighborhood_size = 60,
                     percentile_breaks = c(0,5,10,15,85,90,95,100),
                     neighbor_significance_threshold = 0.9, cluster_min_size = 50,
                     relative_cluster_distance = 30, file_output_prefix = NULL,
                     use_MEM = TRUE, max_alloc = 200000, plot_intensities = FALSE)
{
  require(flowCore)
  require(ggplot2)
  require(ggrastr)
  require(dbscan)
  require(ggrepel)
  require(viridis)
  require(FNN)
  require(grid)
  require(gridExtra)

  # workflow source: https://github.com/cytolab/T-REX

  if(!is.null(file_output_prefix)) {
    file_output_prefix <- paste0(file_output_prefix,"_")
  }
  outdir <- gsub("/$","",outdir)
  if(any(length(compare_list)!=2,class(compare_list)!="list")) {
    stop("error in argument 'compare_list': object must be a named list of length 2. Names of list entries will be used as the categories of comparison. Entries must each be a character vector of file names, with or without file path. File names must match the file names in 'fcs_join_obj' either directly or partially. See unique(fcs_join_obj[['source']])")
  }
  set1_label <- names(compare_list)[1]
  set2_label <- names(compare_list)[2]
  set1_full <- compare_list[[1]]
  set2_full <- compare_list[[2]]
  total_data <- as.data.frame(fcs_join_obj[["data"]])
  total_data$source <- fcs_join_obj[["source"]]
  source_names <- total_data$source
  if(mean(gsub("^.+/","",unlist(compare_list)) %in% unique(source_names))!=0) {
    if(tolower(reduction)=="umap") {
      if("umap" %in% tolower(names(fcs_join_obj))) {
        get_reduction <- fcs_join_obj[["umap"]][["coordinates"]]
      } else {
        stop("error in function call: UMAP specified but no UMAP coordinates found. Run fcs_reduce_dimensions using algorithm = 'UMAP' first")
      }
    } else if(tolower(reduction)=="tsne") {
      if("tsne" %in% tolower(names(fcs_join_obj))) {
        get_reduction <- fcs_join_obj[["tsne"]][["coordinates"]]
      } else {
        stop("error in function call: tSNE specified but no tSNE coordinates found. Run fcs_reduce_dimensions using algorithm = 'tSNE' first")
      }
    }
  } else {
    stop("One or more IDs listed in 'compare_list' were not found in list object's source element.\nHas FCSimple::fcs_reduce_dimensions been run on object?\nIf so, verify that the IDs in 'compare_list' match those in the list object's source element.")
  }
  set1 <- gsub("^.+/","",set1_full)
  set2 <- gsub("^.+/","",set2_full)
  if(length(which(total_data$source %in% gsub("^.+/","",set1)))==0) {
    stop("error in argument 'compare_list': no file names in compare_list[[1]] match data set")
  }
  if(length(which(total_data$source %in% gsub("^.+/","",set2)))==0) {
    stop("error in argument 'compare_list': no file names in compare_list[[2]] match data set")
  }
  set1_ind <- which(total_data$source %in% gsub("^.+/","",set1))
  set2_ind <- which(total_data$source %in% gsub("^.+/","",set2))
  other_ind <- which(!total_data$source %in% gsub("^.+/","",c(set1,set2)))
  set_blank <- rep(NA,times=nrow(total_data))
  set_blank[set1_ind] <- set1_label; set_blank[set2_ind] <- set2_label
  if(length(other_ind)!=0){
    set_blank[other_ind] <- "other"
  }
  total_data$set <- set_blank

  dimred_data <- cbind(total_data, as.data.frame(get_reduction))
  other_ind <- which(!1:nrow(dimred_data) %in% c(set1_ind,set2_ind))

  if(length(other_ind)!=0) {
    set_list <- list(set1 = dimred_data[set1_ind,], set2 = dimred_data[set2_ind,], other = dimred_data[other_ind,])
  } else {
    set_list <- list(set1 = dimred_data[set1_ind,], set2 = dimred_data[set2_ind,])
  }
  for(i in 1:length(set_list)) {
    set_list[[i]] <- split(x = set_list[[i]], f = set_list[[i]]$source)
  }
  set1_size <- sum(sapply(set_list[[1]],nrow)); set2_size <- sum(sapply(set_list[[2]],nrow))
  if(set1_size>set2_size) {
    sample_ratio <- set2_size/set1_size
    for(i in 1:length(set_list[[1]])) {
      set.seed(123)
      set_list[[1]][[i]] <- set_list[[1]][[i]][sample(1:nrow(set_list[[1]][[i]]),size=round(nrow(set_list[[1]][[i]])*sample_ratio),replace=F),]
    }
  } else {
    sample_ratio <- set1_size/set2_size
    for(i in 1:length(set_list[[2]])) {
      set.seed(123)
      set_list[[2]][[i]] <- set_list[[2]][[i]][sample(1:nrow(set_list[[2]][[i]]),size=round(nrow(set_list[[2]][[i]])*sample_ratio),replace=F),]
    }
  }
  ds1 <- do.call(rbind, set_list[[1]]); ds2 <- do.call(rbind, set_list[[2]])
  join_data <- rbind(ds1,ds2)

  if(tolower(reduction)=="umap") {
    map_set <- join_data[,c("UMAP1","UMAP2","set")]
  } else if(tolower(reduction)=="tsne") {
    map_set <- join_data[,c("tSNE1","tSNE2","set")]
  }

  set.seed(123)
  rndmize <- map_set[sample(1:nrow(map_set),size=nrow(map_set),replace=FALSE),]

  pl_pal <- c("blue","red"); names(pl_pal) <- names(compare_list)
  if(tolower(reduction)=="umap") {
    pl <- ggplot(rndmize, aes(x = UMAP1, y = UMAP2, color = set))
  } else if(tolower(reduction)=="tsne") {
    pl <- ggplot(rndmize, aes(x = tSNE1, y = tSNE2, color = set))
  }
  pl <- pl +
    geom_point_rast(pch = 19, size = 0.8, alpha = point_alpha) +
    guides(color = guide_legend(override.aes = list(size = 6, alpha = 1))) +
    scale_color_manual(values = pl_pal) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_text(size = 16),
          legend.position = "bottom")

  ggsave(filename = paste0(file_output_prefix,ifelse(tolower(reduction)=="umap","UMAP_","tSNE_trex_"),strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"), plot = pl,
         device = "pdf", path = outdir, width = 10, height = 10, units = "in", dpi = 900)

  if(tolower(reduction)=="umap") {
    map <- join_data[,c("UMAP1","UMAP2","set")]
  } else if(tolower(reduction)=="tsne") {
    map <- join_data[,c("tSNE1","tSNE2","set")]
  }
  search <- knn.index(data = map[,1:2], k = neighborhood_size)

  search_copy <- search
  search_copy[search_copy %in% 1:nrow(ds1)] <- 0
  search_copy[search_copy!=0] <- 1

  row_mean <- rowMeans(search_copy)
  bin_values <- rep(NA,length(row_mean))

  percentile_breaks <- percentile_breaks[order(percentile_breaks)]
  num_bins <- length(percentile_breaks)-1
  bins <- vector("list", length = num_bins)
  for(i in 1:length(bins)) {
    bins[[i]] <- c(percentile_breaks[i],percentile_breaks[i+1])
    names(bins)[i] <- paste0(bins[[i]],collapse="-")
  }
  bin_index_list <- vector("list", length = length(bins))
  names(bin_index_list) <- names(bins)

  for(i in 1:length(bins)) {
    if(length(bins[[i]])!=2){
      stop("each entry of bins must be of length 2, example: c(1,2)")
    }
    minval <- min(bins[[i]])
    maxval <- max(bins[[i]])
    if(minval==maxval) {
      stop("each entry of bins must consist of two DIFFERENT values that define a range between 0 and 100")
    }
    names(bins)[i] <- paste0(minval,"-",maxval)
    if(maxval==100) {
      bin_index <- intersect(which(row_mean>=minval/100),which(row_mean<=maxval/100))
    } else {
      bin_index <- intersect(which(row_mean>=minval/100),which(row_mean<maxval/100))
    }
    bin_values[bin_index] <- names(bins)[i]
  }

  if(tolower(reduction)=="umap") {
    plot_bins <- map_set[,c("UMAP1","UMAP2","set")]
  } else if(tolower(reduction)=="tsne") {
    plot_bins <- map_set[,c("tSNE1","tSNE2","set")]
  }
  plot_bins$bin <- bin_values
  bin_palette <- colorRampPalette(colors = c("blue","red"))(length(bins))
  if(length(bin_palette)%%2!=0) {
    bin_palette[ceiling(length(bin_palette)/2)] <- "#8C8C8C"
  }
  names(bin_palette) <- names(bins)
  for(i in 1:length(bin_index_list)) {
    bin_index_list[[i]] <- which(plot_bins$bin==names(bin_index_list)[i])
  }

  draw_order <- c()
  pull_from <- 1:length(bin_index_list)
  count <- 1
  repeat {
    if(count%%2==1) {
      sample_index <- pull_from[length(pull_from)]
      pull_from <- pull_from[-length(pull_from)]
    } else {
      sample_index <- pull_from[1]
      pull_from <- pull_from[-1]
    }
    count <- count + 1
    draw_order <- append(draw_order, sample_index)
    if(length(pull_from)==0) {
      break
    }
  }
  point_order <- rev(draw_order)

  draw_index <- c()
  for(i in 1:length(point_order)) {
    draw_index <- append(draw_index, bin_index_list[[point_order[i]]])
  }
  plot_draw_data <- plot_bins[draw_index,]
  if(tolower(reduction)=="umap") {
    pl_bins <- ggplot(data = plot_draw_data, mapping = aes(x = UMAP1, y = UMAP2, color = bin))
  } else if(tolower(reduction)=="tsne") {
    pl_bins <- ggplot(data = plot_draw_data, mapping = aes(x = tSNE1, y = tSNE2, color = bin))
  }
  pl_bins <- pl_bins +
    geom_point_rast(pch = 19, size = 1.5, alpha = point_alpha) +
    scale_color_manual(values = bin_palette, breaks = names(bin_palette)) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(caption = paste0("0 trends toward ",set1_label," and 100 trends toward ",set2_label)) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 12),
          plot.caption = element_text(size = 12, hjust = 0.5))

  ggsave(filename = paste0(file_output_prefix,ifelse(tolower(reduction)=="umap","UMAP","tSNE"),"_trex_by_category_",
                                            strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"), plot = pl_bins,
         device = "pdf", path = outdir, width = 10, height = 10, units = "in", dpi = 900)

  get_high <- which(row_mean>neighbor_significance_threshold)
  get_low <- which(row_mean<(1-neighbor_significance_threshold))

  map_set$bin <- rep(NA,times=nrow(map_set))
  map_set$bin[get_low] <- paste0(set1_label," - ",neighbor_significance_threshold * 100,"%")
  map_set$bin[get_high] <- paste0(set2_label," - ",neighbor_significance_threshold * 100,"%")
  map_set$bin[which(!1:nrow(map_set) %in% c(get_high,get_low))] <- "not significant"
  region_table <- table(map_set$bin)
  write.csv(x = data.frame(regions = region_table),
            file = file.path(outdir,paste0(file_output_prefix,ifelse(tolower(reduction)=="umap","UMAP","tSNE"),"_trex_regions_found_",
                          strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".csv")), row.names = TRUE)
  if(length(region_table)==1) {
    print("no significant regions found. Exiting function.")
    return(0)
  } else {
    print(paste0("regions identified: "))
    print(region_table)
    print(paste0("minimum cluster size was set to: ",cluster_min_size))
    if(max(region_table[-which(names(region_table)=="not significant")])<cluster_min_size) {
      print(paste0("All regions of significance are below min cluster size. Exiting function."))
      return(0)
    }
  }

  map_set_copy <- map_set
  map_set_copy$bin <- factor(x = map_set_copy$bin,
                             levels = c(paste0(set1_label," - ",neighbor_significance_threshold * 100,"%"),
                                        paste0(set2_label," - ",neighbor_significance_threshold * 100,"%"),
                                        "not significant"))

  total_data <- cbind(join_data, map_set_copy)

  map_set_copy <- map_set_copy[order(map_set_copy$bin, decreasing = TRUE),]

  over_cols <- c("blue","red","grey")
  names(over_cols) <- c(paste0(set1_label," - ",neighbor_significance_threshold * 100,"%"),
                        paste0(set2_label," - ",neighbor_significance_threshold * 100,"%"),
                        "not significant")

  if(tolower(reduction)=="umap") {
    pl_data1 <- data.frame(UMAP1 = total_data[,grep("UMAP1",colnames(total_data))[1]],
                           UMAP2 = total_data[,grep("UMAP2",colnames(total_data))[1]],
                           source = total_data$bin)
    pl_hl <- ggplot(data = pl_data1, mapping = aes(x = UMAP1, y = UMAP2, color = source))
  } else if(tolower(reduction)=="tsne") {
    pl_data1 <- data.frame(tSNE1 = total_data[,grep("tSNE1",colnames(total_data))[1]],
                           tSNE2 = total_data[,grep("tSNE2",colnames(total_data))[1]],
                           source = total_data$bin)
    pl_hl <- ggplot(data = pl_data1, mapping = aes(x = tSNE1, y = tSNE2, color = source))
  }
  pl_hl <- pl_hl +
    geom_point_rast(pch = 19, size = 0.7, alpha = 0.5) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    scale_color_manual(values = over_cols) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "bottom")

  ggsave(filename = paste0(file_output_prefix,ifelse(tolower(reduction)=="umap","UMAP","tSNE"),"_trex_significant_",
                           strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"), plot = pl_hl,
         device = "pdf", path = outdir, width = 10, height = 10, units = "in", dpi = 900)

  if(length(which(total_data$bin==paste0(set1_label," - ",neighbor_significance_threshold * 100,"%")))==0) {
    set1_spots <- total_data[1,]
    set1_spots$cluster <- ""
    set1_spots <- set1_spots[-1,]
    set1_dist <- 0
  } else {
    set1_spots <- total_data[which(total_data$bin==paste0(set1_label," - ",neighbor_significance_threshold * 100,"%")),]
    set1_range_x <- diff(range(map_set_copy[,grep("tSNE1|UMAP1",colnames(map_set_copy))[1]]))
    set1_range_y <- diff(range(map_set_copy[,grep("tSNE2|UMAP2",colnames(map_set_copy))[1]]))
    set1_dist <- mean(set1_range_x,set1_range_y)/relative_cluster_distance
  }

  if(length(which(total_data$bin==paste0(set2_label," - ",neighbor_significance_threshold * 100,"%")))==0) {
    set2_spots <- total_data[1,]
    set2_spots$cluster <- ""
    set2_spots <- set2_spots[-1,]
    set2_dist <- 0
  } else {
    set2_spots <- total_data[which(total_data$bin==paste0(set2_label," - ",neighbor_significance_threshold * 100,"%")),]
    set2_range_x <- diff(range(map_set_copy[,grep("tSNE1|UMAP1",colnames(map_set_copy))[1]]))
    set2_range_y <- diff(range(map_set_copy[,grep("tSNE2|UMAP2",colnames(map_set_copy))[1]]))
    set2_dist <- mean(set2_range_x,set2_range_y)/relative_cluster_distance
  }

  set_ns <- total_data[which(total_data$bin=="not significant"),]
  if(nrow(set1_spots)>max_alloc) {
    set.seed(123)
    set1_spots <- set1_spots[sample(1:nrow(set1_spots),size=max_alloc,replace=FALSE),]
    cluster_min_size1 <- cluster_min_size * (max_alloc/nrow(set1_spots))
  } else {
    cluster_min_size1 <- cluster_min_size
  }
  if(nrow(set2_spots)>max_alloc) {
    set.seed(123)
    set2_spots <- set2_spots[sample(1:nrow(set2_spots),size=max_alloc,replace=FALSE),]
    cluster_min_size2 <- cluster_min_size * (max_alloc/nrow(set2_spots))
  } else {
    cluster_min_size2 <- cluster_min_size
  }

  if(set1_dist!=0) {
    if(tolower(reduction)=="umap") {
      scan_set1_cluster <- dbscan(x = set1_spots[,c("UMAP1","UMAP2")],
                                  eps = set1_dist, minPts = cluster_min_size1)$cluster
    } else if(tolower(reduction)=="tsne") {
      scan_set1_cluster <- dbscan(x = set1_spots[,c("tSNE1","tSNE2")],
                                  eps = set1_dist, minPts = cluster_min_size1)$cluster
    }
    set1_spots$cluster <- paste0(set1_label,"_",scan_set1_cluster)
  }

  if(set2_dist!=0) {
    if(tolower(reduction)=="umap") {
      scan_set2_cluster <- dbscan(x = set2_spots[,c("UMAP1","UMAP2")],
                                  eps = set2_dist, minPts = cluster_min_size2)$cluster
    } else if(tolower(reduction)=="tsne") {
      scan_set2_cluster <- dbscan(x = set2_spots[,c("tSNE1","tSNE2")],
                                  eps = set2_dist, minPts = cluster_min_size2)$cluster
    }
    set2_spots$cluster <- paste0(set2_label,"_",scan_set2_cluster)
  }

  set_ns$cluster <- rep("ns_0",nrow(set_ns))

  clustered_data <- do.call(rbind, list(set1_spots, set2_spots, set_ns))
  clustered_data$cluster[grep("0$",clustered_data$cluster)] <- "ns"
  usrc <- unique(clustered_data$source); uclus <- unique(clustered_data$cluster)
  freq_mat <- matrix(data = NA, nrow = length(usrc), ncol = length(uclus))
  row.names(freq_mat) <- usrc; colnames(freq_mat) <- uclus
  for(i in 1:nrow(freq_mat)) {
    tmp_values <- clustered_data$cluster[which(clustered_data$source==row.names(freq_mat)[i])]
    for(j in 1:ncol(freq_mat)) {
      freq_mat[i,j] <- mean(tmp_values==colnames(freq_mat)[j]) * 100
    }
  }
  write.csv(x = freq_mat, file = file.path(outdir,paste0(outdir,"/",file_output_prefix,"trex_significant_cluster_frequencies_",
                                        strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".csv")), row.names = TRUE)

  require(CATALYST)
  require(ComplexHeatmap)
  require(circlize)
  require(grid)
  require(RColorBrewer)

  number_of_cols <- ncol(fcs_join_obj[["data"]])
  event_sources = 1:nrow(clustered_data)
  cluster_numbers = clustered_data$cluster
  include_channels = colnames(clustered_data)[1:number_of_cols]
  heatmap_data = clustered_data[,1:number_of_cols]

  scaled.global <- CATALYST:::.scale_exprs(t(heatmap_data[,include_channels]), 1, 0.01)
  global.t <- t(scaled.global)
  backend.matrix <- matrix(data=NA,nrow=length(unique(cluster_numbers)),ncol=ncol(global.t))
  cluster_numbers <- as.character(cluster_numbers)
  row.names(backend.matrix) <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  colnames(backend.matrix) <- colnames(global.t)
  for(i in 1:nrow(backend.matrix)) {
    get.clus <- which(cluster_numbers==row.names(backend.matrix)[i])
    for(j in 1:ncol(backend.matrix)){
      backend.matrix[i,j] <- median(global.t[get.clus,j])
    }
  }
  hm_pal = rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  z <- backend.matrix
  color.map.fun = circlize::colorRamp2(seq(min(z),max(z), l = n <- 100), colorRampPalette(hm_pal)(n))
  ncell <- rep(NA,times=length(unique(cluster_numbers)))
  names(ncell) <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  for(i in 1:length(ncell)) {
    ncell[i] <- sum(cluster_numbers==names(ncell)[i])
  }
  pop.freq <- matrix(data=ncell,ncol=1)
  row.names(pop.freq) <- names(ncell)
  size_anno_nums <- round((pop.freq/sum(pop.freq))*100,2)
  ranno1 <- rowAnnotation(`Cluster\nSize`=anno_barplot(pop.freq,border=F,width=unit(1.75, "cm"),
                                                       axis_param=list(gp=gpar(fontsize=9)), axis = TRUE),
                          annotation_name_gp=gpar(fontsize=10,fontface="bold"), name = "Cluster\nSize")
  ranno2 <- rowAnnotation(frequency=anno_text(paste0(size_anno_nums,"%"),
                                              gp=gpar(fontsize=10,fontface="bold")))
  backend.matrix <- backend.matrix[order(row.names(backend.matrix)),]
  heatmap_output <- Heatmap(backend.matrix,col=color.map.fun,
                            row_names_side="left",
                            name="median\nscaled\nexpression",
                            heatmap_legend_param=list(at=c(0,0.2,0.4,0.6,0.8,1),legend_height=unit(3,"cm"),
                                                      grid_width=unit(0.6,"cm"),title_position="topleft",
                                                      labels_gp=gpar(fontsize=11),title_gp=gpar(fontsize=11)),
                            row_names_gp=gpar(fontsize=13,fontface="bold"),column_names_gp=gpar(fontsize=12,fontface="bold"),
                            row_gap=unit(1,"mm"),column_gap=unit(1,"mm"),row_dend_gp=gpar(lwd=1.2),row_dend_width=unit(1,"cm"),
                            column_dend_gp = gpar(lwd=1.2), column_dend_height = unit(1,"cm")) +
    ranno1 + ranno2

  ggsave(filename = paste0(file_output_prefix,ifelse(tolower(reduction)=="umap","UMAP","tSNE"),
                           "_trex_significant_cluster_heatmap_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
         plot = grid::grid.grabExpr(draw(heatmap_output)), device = "pdf",
         path = outdir, width = ncol(backend.matrix)/2 + 1,
         height = nrow(backend.matrix)/2 + 1.5, units = "in", dpi = 900, limitsize = FALSE)

  if(use_MEM) {
    # require(cytoMEM)
    require(MEM)
    mem_input <- cbind(heatmap_data,data.frame(cluster = factor(clustered_data$cluster)))
    mem_input$cluster <- as.numeric(mem_input$cluster)
    match_clusters <- data.frame(descriptive_cluster = clustered_data$cluster,
                                 numeric_cluster = mem_input$cluster)
    match_clusters <- match_clusters[-which(duplicated(match_clusters$descriptive_cluster)),]
    clus_order <- order(match_clusters$numeric_cluster)

    mcalc <- MEM::MEM(exp_data = mem_input, transform=FALSE, choose.markers=FALSE,
                      rename.markers=FALSE, choose.ref=FALSE)

    check_mem <- sapply(mcalc,function(x) return(nrow(x[[1]])))
    for(i in 1:length(check_mem)) {
      if(is.null(unlist(check_mem[i]))) {
        next
      }
      if(check_mem[i]==nrow(match_clusters)) {
        mcalc[[i]][[1]] <- mcalc[[i]][[1]][clus_order,]
        row.names(mcalc[[i]][[1]]) <- match_clusters$descriptive_cluster[clus_order]
      }
    }

    MEM::build.heatmaps(mcalc, cluster.MEM = "none", cluster.medians = "none",
                        display.thresh = 1,  output.files = TRUE, labels = FALSE,
                        only.MEMheatmap = TRUE)
  }
  mem_outs <- list.files(path = file.path(outdir,paste0(getwd(),"/output files")), full.names = TRUE)
  for(i in 1:length(mem_outs)) {
    file.copy(from = mem_outs[i], to = outdir, overwrite = TRUE,
              recursive = FALSE, copy.mode = TRUE)
    file.remove(mem_outs[i])
  }

  if(tolower(reduction)=="umap") {
    plot_data <- clustered_data[,c("UMAP1","UMAP2","cluster")]
  } else if(tolower(reduction)=="tsne") {
    plot_data <- clustered_data[,c("tSNE1","tSNE2","cluster")]
  }
  sig_clus <- plot_data[-which(plot_data$cluster=="ns"),]
  xclus <- rep(NA,length(unique(sig_clus$cluster))); names(xclus) <- unique(sig_clus$cluster); yclus <- xclus
  if(tolower(reduction)=="umap"){
    for(i in 1:length(xclus)) {
      xclus[i] <- median(sig_clus$UMAP1[which(sig_clus$cluster==names(xclus)[i])])
      yclus[i] <- median(sig_clus$UMAP2[which(sig_clus$cluster==names(yclus)[i])])
    }
  } else if(tolower(reduction)=="tsne"){
    for(i in 1:length(xclus)) {
      xclus[i] <- median(sig_clus$tSNE1[which(sig_clus$cluster==names(xclus)[i])])
      yclus[i] <- median(sig_clus$tSNE2[which(sig_clus$cluster==names(yclus)[i])])
    }
  }
  if(tolower(reduction)=="umap") {
    pl_lab <- ggplot(data = plot_data, mapping = aes(x = UMAP1, y = UMAP2, color = cluster))
  } else if(tolower(reduction)=="tsne") {
    pl_lab <- ggplot(data = plot_data, mapping = aes(x = tSNE1, y = tSNE2, color = cluster))
  }
  pl_lab <- pl_lab +
    geom_point_rast(pch = 19, size = 0.5, alpha = 0.1) +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    annotate("text_repel", x = xclus, y = yclus, label = names(xclus), size = 5) +
    theme_void() +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "none")

  sig_only <- plot_data[-which(plot_data$cluster=="ns"),]
  if(tolower(reduction)=="umap") {
    pl_sig_lab <- ggplot(data = sig_only, mapping = aes(x = UMAP1, y = UMAP2, color = cluster)) +
      geom_point_rast(pch = 19, size = 0.7, alpha = 0.5) +
      xlim(range(plot_data$UMAP1)) + ylim(range(plot_data$UMAP2))
  } else if(tolower(reduction)=="tsne") {
    pl_sig_lab <- ggplot(data = sig_only, mapping = aes(x = tSNE1, y = tSNE2, color = cluster)) +
      geom_point_rast(pch = 19, size = 0.7, alpha = 0.5) +
      xlim(range(plot_data$tSNE1)) + ylim(range(plot_data$tSNE2))
  }
  pl_sig_lab <- pl_sig_lab +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    annotate("text_repel", x = xclus, y = yclus, label = names(xclus), size = 5) +
    theme_void() +
    theme(legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          legend.text = element_text(size = 14),
          legend.position = "none")

  listed_plots <- mget(c("pl_lab","pl_sig_lab"))
  ggsave(filename = paste0(ifelse(tolower(reduction)=="umap","UMAP","tSNE"),"_trex_significant_labeled_",
                           strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
         plot = gridExtra::arrangeGrob(grobs = listed_plots, nrow=2, ncol=1),
         device = "pdf", path = outdir, width = 10, height = 20, units = "in", dpi = 900)

  if(plot_intensities) {
    intens_list <- vector("list", length = number_of_cols); names(intens_list) <- colnames(join_data)[1:number_of_cols]
    for(i in 1:length(intens_list)) {
      intens_list[[i]] <- cbind(join_data[,names(intens_list)[i]],map); colnames(intens_list[[i]])[1] <- names(intens_list)[i]
    }
    intens_pl <- function(arg1, method = "color") {
      if(nrow(arg1)>200000){
        arg1 <- arg1[sample(1:nrow(arg1),size=200000,replace=F),]
      }
      cname1 <- colnames(arg1)[1]
      colnames(arg1)[1] <- "col1"
      if(method=="color") {
        if(tolower(reduction)=="umap"){
          plt <- ggplot(data = arg1, mapping = aes_string(x = "UMAP1", y = "UMAP2", color = colnames(arg1)[1]))
        } else if(tolower(reduction)=="tsne") {
          plt <- ggplot(data = arg1, mapping = aes_string(x = "tSNE1", y = "tSNE2", color = colnames(arg1)[1]))
        }
        plt <- plt +
          geom_point_rast(size = 0.8, pch = 19, alpha = 0.1) +
          scale_color_viridis(option = "D") +
          ggtitle(cname1) +
          theme_minimal() +
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                legend.title = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.position = "bottom")
      } else if(method=="fill") {
        if(tolower(reduction)=="umap") {
          plt <- ggplot(data = arg1, mapping = aes_string(x = "UMAP1", y = "UMAP2", fill = colnames(arg1)[1]))
        } else if(tolower(reduction)=="tsne") {
          plt <- ggplot(data = arg1, mapping = aes_string(x = "tSNE1", y = "tSNE2", fill = colnames(arg1)[1]))
        }
        plt <- plt +
          geom_point_rast(size = 0.8, pch = 21, stroke = 0.05) +
          scale_fill_viridis(option = "D") +
          ggtitle(cname1) +
          theme_minimal() +
          theme(axis.title = element_blank(),
                axis.text = element_blank(),
                legend.title = element_blank(),
                plot.title = element_text(hjust = 0.5),
                legend.position = "bottom")
      }
      return(plt)
    }
    intens_plots <- lapply(intens_list, intens_pl)

    void_plot <- ggplot(data = cars, aes(x=speed,y=dist)) +
      geom_point(alpha = 0) +
      theme_void()

    arrange_length <- ceiling(length(intens_list)/4)
    for(i in 1:(arrange_length*4)) {
      if(i>length(intens_plots)) {
        intens_plots[[i]] <- void_plot
      }
    }

    arranged_list <- vector("list", length = arrange_length)
    low_index <- seq(from = 1, to = arrange_length*4, by = 4)
    high_index <- seq(from = 4, to = arrange_length*4, by = 4)

    for(i in 1:length(arranged_list)) {
      arranged_list[[i]] <- ggpubr::ggarrange(plotlist = intens_plots[low_index[i]:high_index[i]],
                                              nrow = 2, ncol = 2)
    }

    ggsave(filename = paste0(file_output_prefix,ifelse(tolower(reduction)=="umap","UMAP","tSNE"),"_trex_heatmaps_",
                             strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
           plot = gridExtra::marrangeGrob(grobs = arranged_list, nrow=1, ncol=1, top = ""),
           device = "pdf", path = outdir, width = 12, height = 12, units = "in", dpi = 900)
  }
  # return(fcs_join_obj)
}
