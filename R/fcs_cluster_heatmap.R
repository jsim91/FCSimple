fcs_cluster_heatmap <- function(fcs_join_obj, algorithm = "leiden",
                                heatmap_color_palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                fold_method = "median", heatmap_intensity_method = "median",
                                tile_data = NULL, include_colnames = TRUE,
                                cl_rows = TRUE, cl_cols = TRUE, show_legend = TRUE,
                                scale_local = NULL,
                                fold_method_name = "difference", group_names = NULL,
                                heatmap_title = "")
{
  require(CATALYST)
  require(ComplexHeatmap)
  require(circlize)

  scaled.global <- CATALYST:::.scale_exprs(t(heatmap_data[,include_channels]), 1, 0.01)
  global.t <- t(scaled.global)
  backend.matrix <- matrix(data=NA,nrow=length(unique(cluster_numbers)),ncol=ncol(global.t))
  cluster_numbers <- as.character(cluster_numbers)
  row.names(backend.matrix) <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  colnames(backend.matrix) <- colnames(global.t)
  for(i in 1:nrow(backend.matrix)){
    get.clus <- which(cluster_numbers==row.names(backend.matrix)[i])
    for(j in 1:ncol(backend.matrix)){
      backend.matrix[i,j] <- median(global.t[get.clus,j])
    }
  }
  hm_pal = heatmap_color_palette
  z <- backend.matrix
  color.map.fun = circlize::colorRamp2(seq(min(z),max(z), l = n <- 100), colorRampPalette(hm_pal)(n))
  ncell <- rep(NA,times=length(unique(cluster_numbers)))
  names(ncell) <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  for(i in 1:length(ncell)){
    ncell[i] <- sum(cluster_numbers==names(ncell)[i])
  }
  pop.freq <- matrix(data=ncell,ncol=1)
  row.names(pop.freq) <- names(ncell)
  size_anno_nums <- round((pop.freq/sum(pop.freq))*100,2)
  # ranno1 <- rowAnnotation(`Cluster\nSize`=anno_barplot(pop.freq,border=F,width=unit(1.75, "cm"),
  #                                                      axis_param=list(gp=gpar(fontsize=9)), axis = TRUE),
  #                         annotation_name_gp=gpar(fontsize=10,fontface="bold"), name = "Cluster\nSize")
  ranno1 <- rowAnnotation(`Cluster\nSize`=anno_barplot(pop.freq,border=F,width=unit(1.75, "cm"),
                                                       axis_param=list(gp=gpar(fontsize=9), at = 0, labels = ""),
                                                       axis = TRUE, just = "center", location = unit(0.5, "npc"), show_name = TRUE),
                          annotation_name_gp=gpar(fontsize=10,fontface="bold"), annotation_name_rot = 0, name = "Cluster\nSize")
  ranno2 <- rowAnnotation(frequency=anno_text(paste0(size_anno_nums,"%"),
                                              gp=gpar(fontsize=10,fontface="bold")))
  heatmap_output <- Heatmap(backend.matrix,col=color.map.fun,
                            row_names_side="left",
                            name="median\nscaled\nexpression",
                            heatmap_legend_param=list(at=c(0,0.2,0.4,0.6,0.8,1),legend_height=unit(3,"cm"),
                                                      grid_width=unit(0.6,"cm"),title_position="topleft",
                                                      labels_gp=gpar(fontsize=11),title_gp=gpar(fontsize=11,fontface="bold")),
                            row_names_gp=gpar(fontsize=13,fontface="bold"),column_names_gp=gpar(fontsize=12,fontface="bold"),
                            row_gap=unit(1,"mm"),column_gap=unit(1,"mm"),row_dend_gp=gpar(lwd=1.5),row_dend_width=unit(1,"cm"),
                            column_dend_gp = gpar(lwd=1.5), column_dend_height = unit(1,"cm")) +
    ranno1 + ranno2
  return(list(heatmap = heatmap_output, heatmap_tile_data = backend.matrix, population_size = pop.freq))
}
