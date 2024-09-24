fcs_plot_reduction_difference <- function(fcs_join_obj, reduction = c("UMAP","tSNE"),
                                          compare_list, color_list, n_kde = 200,
                                          outdir = getwd(), axis_title_text_size = 12,
                                          dbscan_eps = "auto", dbscan_minPts = "auto",
                                          legend_label_text_size = 12, annotate_clusters = TRUE,
                                          cluster_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                          cluster_number_annotation_size = 5, contour_bin_width = 0.001,
                                          legend_orientation = c("horizontal","vertical"),
                                          figure_width = 8, figure_height = 8, hull_radius = 1.5,
                                          add_timestamp = TRUE, hull_concavity = 2)
{
  require(ggplot2)
  require(MASS)
  require(reshape2)
  require(scales)
  require(ggrastr)
  require(ggforce); require(concaveman)
  require(dbscan)
  require(ggpubr)
  require(shadowtext)
  require(ComplexHeatmap)

  # testing
  # fcs_join_obj = fcs_obj
  # reduction = "UMAP"
  # compare_list = compare_hiv
  # color_list = color_hiv
  # n_kde = 200
  # outdir = paste0(getwd(),"/HIV_status_test")
  # axis_title_text_size = 12
  # hull_radius = 2
  # dbscan_eps = "auto"
  # dbscan_minPts = "auto"
  # legend_label_text_size = 12
  # annotate_clusters = TRUE
  # cluster_algorithm = "flowsom"
  # cluster_number_annotation_size = 5
  # contour_bin_width = 0.001
  # legend_orientation = "vertical"
  # figure_width = 8
  # figure_height = 8
  # add_timestamp = TRUE
  # hull_concavity = 2

  if(length(reduction)!=1) {
    stop("error in argument 'reduction': Use either tSNE or UMAP. Reduction must be present in 'fcs_join_obj'.")
  }
  if(length(legend_orientation)!=1) {
    warning("Use either 'horizontal' or 'vertical'. Proceeding with 'legend_orientation' = 'horizontal'")
    legend_orientation = "horizontal"
  }

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  reduction_names <- colnames(reduction_coords)[1:2]
  colnames(reduction_coords)[1:2] <- c("dimx","dimy")
  dim1_range <- range(reduction_coords[,1]); dim2_range <- range(reduction_coords[,2])

  if(length(cluster_algorithm)!=1) {
    cluster_algorithm <- cluster_algorithm[1]
  }
  if(annotate_clusters) {
    obj_clusters <- as.character(fcs_join_obj[[tolower(cluster_algorithm)]][["clusters"]])
    unique_clus <- unique(obj_clusters); unique_clus <- unique_clus[order(as.numeric(unique_clus))]
    clusx <- rep(NA,length(unique_clus)); names(clusx) <- unique_clus; clusy <- clusx
    for(i in 1:length(clusx)) {
      clusx[i] <- median(reduction_coords[,1][which(obj_clusters==names(clusx)[i])])
      clusy[i] <- median(reduction_coords[,2][which(obj_clusters==names(clusx)[i])])
    }
  }

  if(any(length(compare_list)!=2, length(color_list)!=2)) {
    stop("error in argument 'compare_list' and 'color_list': both lists need to be length 2")
  }
  if(mean(names(compare_list) %in% names(color_list))!=1) {
    stop("error in argument 'compare_list' and 'color_list': list names must be the same in each")
  }
  if(any(length(unique(names(compare_list)))!=2, length(unique(names(color_list)))!=2)) {
    stop("error in argument 'compare_list' and 'color_list': names within each list must be unique")
  }
  grp1_reduction_loc <- which(fcs_join_obj[["source"]] %in% compare_list[[1]])
  grp2_reduction_loc <- which(fcs_join_obj[["source"]] %in% compare_list[[2]])
  grp1_red <- reduction_coords[grp1_reduction_loc,]#; grp1_red <- rbind(grp1_red, pad_mat)
  grp2_red <- reduction_coords[grp2_reduction_loc,]#; grp2_red <- rbind(grp2_red, pad_mat)
  ds_to <- min(nrow(grp1_red),nrow(grp2_red))
  if(nrow(grp1_red)>ds_to) {
    set.seed(123)
    grp1_red <- grp1_red[sample(x = 1:nrow(grp1_red), size = ds_to, replace = FALSE),]
  }
  if(nrow(grp2_red)>ds_to) {
    set.seed(123)
    grp2_red <- grp2_red[sample(x = 1:nrow(grp2_red), size = ds_to, replace = FALSE),]
  }

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
  colnames(diff12.m)[3] <- "z"

  minlim <- min(diff12.m$z); maxlim <- max(diff12.m$z)

  # where group2 is list element 2, positive values indicate some enrichment for group2

  set.seed(123)
  background_data <- reduction_coords[sample(1:nrow(reduction_coords),size=ifelse(nrow(reduction_coords)>100000,100000,nrow(reduction_coords)),replace=FALSE),]

  plt_dens_diff <- ggplot(diff12.m, aes(x = Var1, y = Var2, z=z, fill=z)) +
    geom_tile() +
    # stat_contour(aes(colour = after_stat(!!str2lang("level"))), binwidth = contour_bin_width) +
    stat_contour(aes(colour = after_stat(level)), binwidth = contour_bin_width) +
    scale_fill_gradient2(low = color_list[[1]],mid = "white",
                         high = color_list[[2]], midpoint = 0, breaks = c(minlim,maxlim),
                         labels=c(names(compare_list)[1],names(compare_list)[2]),limits=c(minlim,maxlim)) +
    scale_colour_gradient2(low = muted(color_list[[1]]), mid = "white", high = muted(color_list[[2]]), midpoint = 0) +
    coord_cartesian(xlim = dim1_range, ylim = dim2_range, expand = FALSE) +
    guides(color = "none", fill = guide_colorbar(ticks.color = NA)) +
    labs(x = reduction_names[1], y = reduction_names[2])
  if(annotate_clusters) {
    plt_dens_diff <- plt_dens_diff +
      annotate("shadowtext", x = clusx, y = clusy, label = names(clusx), size = cluster_number_annotation_size) +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.direction= legend_orientation,
            legend.title = element_blank(),
            legend.text = element_text(size = legend_label_text_size, hjust = 0.5, vjust = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(face = "bold", size = axis_title_text_size))
  } else {
    plt_dens_diff <- plt_dens_diff +
      theme_bw() +
      theme(legend.position = "bottom",
            legend.direction= legend_orientation,
            legend.title = element_blank(),
            legend.text = element_text(size = legend_label_text_size, hjust = 0.5, vjust = 1),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(face = "bold", size = axis_title_text_size))
  }

  col_fun = circlize::colorRamp2(seq(min(diff12.m$z),max(diff12.m$z), l = n <- 100), colorRampPalette(unlist(color_list))(n))
  lgd_h <- Legend(col_fun = col_fun, title = "", border = "black", direction = "horizontal",
                  at = c(min(diff12.m$z), 0, max(diff12.m$z)), legend_width = unit(6, "cm"), legend_height = unit(1, "cm"),
                  title_position = "topcenter", labels = c(names(color_list)[1], "", names(color_list)[2]),
                  labels_gp = gpar(fontsize = 20), tick_length = unit(0,"npc"))
  lgd_v <- Legend(col_fun = col_fun, title = "", border = "black", direction = "vertical",
                  at = c(min(diff12.m$z), 0, max(diff12.m$z)), legend_width = unit(1, "cm"), legend_height = unit(6, "cm"),
                  title_position = "topcenter", labels = c(names(color_list)[1], "", names(color_list)[2]),
                  labels_gp = gpar(fontsize = 20), tick_length = unit(0,"npc"))
  # dev.off() # testing
  # pushViewport(viewport(width = 0.9, height = 0.9))
  # draw(lgd_h, x = unit(0.25, "npc"), y = unit(0.5, "npc"))
  # draw(lgd_v, x = unit(0.7, "npc"), y = unit(0.5, "npc"))
  # popViewport()

  # plt_leg <- ggpubr::as_ggplot(ggpubr::get_legend(p = plt_dens_diff))
  lgh_gg <- as_ggplot(grid.grabExpr(expr = draw(lgd_h, x = unit(0.5, "npc"), y = unit(0.5, "npc"))))
  lgv_gg <- as_ggplot(grid.grabExpr(expr = draw(lgd_v, x = unit(0.5, "npc"), y = unit(0.5, "npc"))))
  lg_arr <- ggpubr::ggarrange(plotlist = list(lgh_gg, lgv_gg), nrow = 1)

  if(length(dbscan_eps)!=1) {
    dbscan_eps <- "auto"
  }
  if(length(dbscan_minPts)!=1) {
    dbscan_minPts <- "auto"
  }
  if(dbscan_eps == "auto") {
    dbscan_eps_arg <- mean(diff(dim1_range),diff(dim2_range))/20
  } else {
    dbscan_eps_arg <- dbscan_eps
  }
  if(dbscan_minPts == "auto") {
    dbscan_minpts_arg <- floor(nrow(background_data)/100)
  } else {
    dbscan_minpts_arg <- dbscan_minPts
  }
  scan_islands <- dbscan::dbscan(x = background_data[,c(1:2)], eps = dbscan_eps_arg, minPts = dbscan_minpts_arg)
  background_data$island <- as.character(scan_islands$cluster)
  background_data$island_1 <- scan_islands$cluster
  override_island_col <- rep(alpha(c("black"), alpha = 1), length(unique(background_data$island)))
  names(override_island_col) <- unique(background_data$island)

  plt_dens_back <- ggplot(background_data, aes(x = dimx, y = dimy, color = island_1)) +
    geom_point(alpha = 0) +
    geom_mark_hull(concavity = hull_concavity, expand = unit(0.25, "mm"),
                   aes(fill = island, filter = island != '0'),
                   color = alpha(colour = "white", alpha = 0), radius = unit(hull_radius, "mm")) +
    scale_fill_manual(values = override_island_col) +
    coord_cartesian(xlim = dim1_range, ylim = dim2_range, expand = FALSE) +
    guides(fill = "none", color = guide_colorbar(ticks.color = NA)) +
    labs(x = reduction_names[1], y = reduction_names[2]) +
    theme_bw() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = legend_label_text_size, hjust = 0.5, vjust = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(face = "bold", size = axis_title_text_size, color = "white"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  if(annotate_clusters) {
    plt_dens_back <- plt_dens_back + annotate("shadowtext", x = clusx, y = clusy, label = names(clusx), size = cluster_number_annotation_size)
  }

  timestamp <- strftime(Sys.time(),"%Y-%m-%d_%H%M%S")
  if(add_timestamp) {
    fname_top <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                        names(compare_list)[2],"_",tolower(reduction),
                        "_density_difference_over_",
                        timestamp)
    fname_bottom <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_",tolower(reduction),
                           "_density_difference_under_",
                           timestamp)
    fname_legend <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_",tolower(reduction),
                           "_density_difference_legend_",
                           timestamp)
    fname <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                    names(compare_list)[2],"_",tolower(reduction),
                    "_density_difference_",
                    timestamp)
  } else {
    fname_top <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                        names(compare_list)[2],"_",tolower(reduction),
                        "_density_difference_over")
    fname_bottom <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_",tolower(reduction),
                           "_density_difference_under")
    fname_legend <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_",tolower(reduction),
                           "_density_difference_legend")
    fname <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                    names(compare_list)[2],"_",tolower(reduction),
                    "_density_difference")
  }

  rm_box <- theme(plot.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank())

  umap_arr <- ggpubr::ggarrange(plotlist = list(plt_dens_back + rm_box, plt_dens_diff + rm_box + theme(legend.position = "none", axis.title = element_blank())), nrow = 1)

  ggsave(filename = paste0(fname_bottom,".pdf"),
         plot = plt_dens_back, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 900, bg = "white")
  ggsave(filename = paste0(fname_top,".pdf"),
         plot = plt_dens_diff + theme(legend.position = "none"), device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 900, bg = "white")
  ggsave(filename = paste0(fname_legend,".pdf"),
         plot = lg_arr, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 900, bg = "white")
  ggsave(filename = paste0(fname,".pdf"),
         plot = umap_arr, device = "pdf",
         width = figure_width*2, height = figure_height,
         units = "in", dpi = 900, bg = "white")
}
