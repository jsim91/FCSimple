fcs_cluster_heatmap <- function(fcs_join_obj, algorithm,
                                heatmap_color_palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
{
  if(!tolower(algorithm) %in% names(fcs_join_obj)) {
    stop("error in argument 'algorithm': algorithm not found in fcs_join_obj. Try 'View(fcs_join_obj)'")
  }

  require(CATALYST)
  require(ComplexHeatmap)
  require(circlize)
  require(grid)

  event_source <- fcs_join_obj[["source"]]
  cluster_numbers <- as.numeric(fcs_join_obj[[tolower(algorithm)]][["clusters"]])
  heatmap_data <- fcs_join_obj[["data"]]
  include_channels <- colnames(fcs_join_obj[["data"]])

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
  hm_pal = heatmap_color_palette
  z <- backend.matrix
  color.map.fun = circlize::colorRamp2(seq(min(z),max(z), l = n <- 100), colorRampPalette(hm_pal)(n))
  ncell <- rep(NA,times=length(unique(cluster_numbers)))
  names(ncell) <- unique(cluster_numbers)[order(as.numeric(unique(cluster_numbers)))]
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
  backend.matrix <- backend.matrix[order(as.numeric(row.names(backend.matrix))),]
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
  fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]] <- list(heatmap = heatmap_output,
                                                                heatmap_tile_data = backend.matrix,
                                                                population_size = pop.freq)
  return(fcs_join_obj)
}

fcs_plot_heatmap <- function(fcs_join_obj, algorithm, outdir = getwd())
{
  require(ggplot2)
  ggsave(filename = paste0(outdir,"/",tolower(algorithm),"_cluster_heatmap_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
         plot = grid::grid.grabExpr(draw(fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap"]])),
         device = "pdf", width = (ncol(fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]])*0.33)+3.25,
         height = (nrow(fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]])*0.3)+2.25,
         units = "in", dpi = 900)
}
