fcs_plot_reduction_difference <- function(fcs_join_obj, reduction = c("UMAP","tSNE"),
                                          compare_list, color_list, n_kde = 200,
                                          outdir = getwd(), axis_title_text_size = 12,
                                          legend_label_text_size = 12)
{
  require(ggplot2)
  require(MASS)
  require(reshape2)
  require(scales)

  # testing
  # fcs_join_obj <- readRDS("J:/umap_cluster_object.rds")
  # reduction <- "umap"
  # grp1_ind <- grep(pattern = "V1", x = fcs_join_obj$source)
  # grp2_ind <- grep(pattern = "V2", x = fcs_join_obj$source)
  # compare_list <- list(visit1 = grp1_ind, visit2 = grp2_ind)
  # color_list <- list(visit1 = "#1bbc9b", visit2 = "#ff822e")
  # n_kde <- 200
  # axis_title_text_size = 12
  # legend_label_text_size = 12
  # outdir = getwd()

  if(length(reduction)!=1) {
    stop("error in argument 'reduction': Use either tSNE or UMAP. Reduction must be present in 'fcs_join_obj'.")
  }

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  dim1_range <- range(reduction_coords[,1]); dim2_range <- range(reduction_coords[,2])
  pad1 <- c(min(dim1_range*1.1),min(dim2_range*1.1)); pad2 <- c(min(dim1_range*1.1),max(dim2_range*1.1))
  pad3 <- c(max(dim1_range*1.1),min(dim2_range*1.1)); pad4 <- c(max(dim1_range*1.1),max(dim2_range*1.1))
  pad_mat <- rbind(rbind(rbind(matrix(data = pad1, nrow = 1, ncol = 2), matrix(data = pad2, nrow = 1, ncol = 2)),
                         matrix(data = pad3, nrow = 1, ncol = 2)),
                   matrix(data = pad4, nrow = 1, ncol = 2))
  colnames(pad_mat) <- colnames(reduction_coords)

  if(any(length(compare_list)!=2, length(color_list)!=2)) {
    stop("error in argument 'compare_list' and 'color_list': both lists need to be length 2")
  }
  if(mean(names(compare_list) %in% names(color_list))!=1) {
    stop("error in argument 'compare_list' and 'color_list': list names must be the same in each")
  }
  if(any(length(unique(names(compare_list)))!=2, length(unique(names(color_list)))!=2)) {
    stop("error in argument 'compare_list' and 'color_list': names within each list must be unique")
  }
  grp1_red <- reduction_coords[compare_list[[1]],]; grp1_red <- rbind(grp1_red, pad_mat)
  grp2_red <- reduction_coords[compare_list[[2]],]; grp2_red <- rbind(grp2_red, pad_mat)

  # method source: https://stackoverflow.com/questions/28521145/r-calculate-and-plot-difference-between-two-density-countours

  grp1_kde2d <- MASS::kde2d(x = grp1_red[,1], y = grp1_red[,2], lims = c(dim1_range, dim2_range), n = n_kde)
  grp2_kde2d <- MASS::kde2d(x = grp2_red[,1], y = grp2_red[,2], lims = c(dim1_range, dim2_range), n = n_kde)
  if(any(!identical(grp1_kde2d$x, grp2_kde2d$x),!identical(grp1_kde2d$y, grp2_kde2d$y))) {
    stop("kde2d x,y values are not identical")
  }
  diff12 = grp1_kde2d
  diff12$z = grp2_kde2d$z - grp1_kde2d$z
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y
  diff12.m = melt(diff12$z, id.var=rownames(diff12))
  names(diff12.m) = c("UMAP1","UMAP2","z")

  minlim <- min(diff12.m$z); maxlim <- max(diff12.m$z)

  # where group2 is list element 2, positive values indicate some enrichment for group2

  plt_dens_diff <- ggplot(diff12.m, aes(UMAP1, UMAP2, z=z, fill=z)) +
    geom_tile() +
    stat_contour(aes(colour=..level..), binwidth = 0.001) +
    scale_fill_gradient2(low = "red",mid = "white", high = "blue", midpoint = 0, breaks = c(minlim,maxlim),
                         labels=c(names(compare_list[1]),names(compare_list[2])),limits=c(minlim,maxlim)) +
    scale_colour_gradient2(low = muted("red"), mid = "white", high = muted("blue"), midpoint = 0) +
    coord_cartesian(xlim = dim1_range, ylim = dim2_range, expand = FALSE) +
    guides(color = "none", fill = guide_colorbar(ticks.color = NA)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = legend_label_text_size, hjust = 0.5, vjust = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(face = "bold", size = axis_title_text_size))

  fname <- paste0(outdir,"/",names(compare_list)[1],"_vs_",names(compare_list)[2],"_reduction_density_difference_",
                  strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))

  ggsave(filename = paste0(fname,".pdf"),
         plot = plt_dens_diff, device = "pdf", width = 8, height = 8,
         units = "in", dpi = 900, bg = "white")
}
