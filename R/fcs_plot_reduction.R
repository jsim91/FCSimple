fcs_plot_reduction <- function(fcs_join_obj, algorithm, reduction, point_alpha = 0.1, outdir = getwd(),
                               internal_call = FALSE, anno_indices = NULL, keep_indices = NA, pdf_dim = 10,
                               png_dim = 1000, plotting_device = "pdf", return_plot = TRUE)
{
  require(ggplot2)
  require(shadowtext)
  require(ggrastr)

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  cluster_numbers <- as.numeric(as.character(fcs_join_obj[[tolower(algorithm)]][["clusters"]]))
  uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  xval <- rep(NA,times=length(uclus)); names(xval) <- uclus; yval <- xval
  for(i in 1:length(xval)) {
    xval[i] <- median(reduction_coords[,1][which(cluster_numbers==as.numeric(names(xval)[i]))])
    yval[i] <- median(reduction_coords[,2][which(cluster_numbers==as.numeric(names(yval)[i]))])
  }
  plt_input <- cbind(reduction_coords,data.frame(cluster = cluster_numbers))
  plt_input$cluster <- factor(plt_input$cluster)
  if(!internal_call) {
    plt_reduction <- ggplot(data = plt_input, mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                   y = colnames(reduction_coords)[2],
                                                                   color = "cluster")) +
      ggrastr::geom_point_rast(alpha = point_alpha) +
      annotate("shadowtext", x = xval, y = yval, label = names(xval), size = 5) +
      theme_void() +
      theme(legend.position = "none")
    fname <- paste0(outdir,"/",tolower(algorithm),"_",tolower(reduction),"_labeled_",
                    strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf")
  } else {
    plt_reduction <- ggplot(data = plt_input[-keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                       y = colnames(reduction_coords)[2])) +
      ggrastr::geom_point_rast(alpha = point_alpha, color = "grey") +
      ggrastr::geom_point_rast(data = plt_input[keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                        y = colnames(reduction_coords)[2]),
                               alpha = point_alpha, color = "red") +
      annotate("shadowtext", x = xval, y = yval, label = names(xval), size = 5) +
      theme_void() +
      theme(legend.position = "none")
    fname <- paste0(outdir,"/islands_selected_for_by_dbscan_",
                    strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
  }
  if(return_plot) {
    return(plt_reduction)
  } else {
    if(plotting_device=="pdf") {
      ggsave(filename = paste0(fname,".pdf"),
             plot = plt_reduction, device = "pdf", width = pdf_dim, height = pdf_dim,
             units = "in", dpi = 900)
    } else if(plotting_device=="png"){
      ggsave(filename = paste0(fname,".png"),
             plot = plt_reduction, device = "png", width = png_dim, height = png_dim,
             units = "in", dpi = 900)
    }
  }
}
