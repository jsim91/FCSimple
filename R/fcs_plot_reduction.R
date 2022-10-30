fcs_plot_reduction <- function(fcs_join_obj, algorithm, reduction, point_alpha = 0.1, outdir = getwd())
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
  plt_reduction <- ggplot(data = plt_input, mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                 y = colnames(reduction_coords)[2],
                                                                 color = "cluster")) +
    ggrastr::geom_point_rast(alpha = point_alpha) +
    annotate("shadowtext", x = xval, y = yval, label = names(xval), size = 5) +
    theme_void()
  ggsave(filename = paste0(outdir,"/",tolower(algorithm),"_",tolower(reduction),"_labeled_",
                           strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
         plot = plt_reduction, device = "pdf", width = 10, height = 10,
         units = "in", dpi = 900)
}
