fcs_plot_reduction_difference <- function(fcs_join_obj, reduction = c("UMAP","tSNE"),
                                          compare_list, color_list, n_kde = 200,
                                          outdir = getwd(), axis_title_text_size = 12,
                                          dbscan_eps_ratio = "auto", dbscan_minPts_ratio = "auto",
                                          legend_label_text_size = 12, annotate_clusters = TRUE,
                                          cluster_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                          cluster_number_annotation_size = 5,
                                          legend_orientation = c("horizontal","vertical"),
                                          figure_width = 8, figure_height = 8,
                                          add_timestamp = TRUE)
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

  # testing
  # fcs_join_obj <- readRDS("J:/Mashayekhi/flow/project_1/outs/fcsimple_join_obj_clustered.rds")
  # reduction <- "umap"
  # grp1_ind <- grep(pattern = "SD1", x = fcs_join_obj$source)
  # grp2_ind <- grep(pattern = "SD2", x = fcs_join_obj$source)
  # compare_list <- list(SD1 = grp1_ind, SD2 = grp2_ind)
  # color_list <- list(SD1 = "#1bbc9b", SD2 = "#ff822e")
  # n_kde <- 200
  # axis_title_text_size = 12
  # legend_label_text_size = 12
  # outdir = getwd()
  # legend_pos = c(0.5,0.075)
  # legend_orientation = "horizontal"
  # add_timestamp = TRUE
  # annotate_clusters = TRUE
  # cluster_algorithm = "leiden"
  # dbscan_eps_ratio = "auto"
  # dbscan_minPts_ratio = "auto"
  # cluster_number_annotation_size = 6

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
    obj_clusters <- fcs_join_obj[[tolower(cluster_algorithm)]][["clusters"]]
    unique_clus <- unique(obj_clusters); unique_clus <- unique_clus[order(as.numeric(unique_clus))]
    clusx <- rep(NA,length(unique_clus)); names(clusx) <- unique_clus; clusy <- clusx
    for(i in 1:length(clusx)) {
      clusx[i] <- median(reduction_coords[,1][which(obj_clusters==as.numeric(names(clusx)[i]))])
      clusy[i] <- median(reduction_coords[,2][which(obj_clusters==as.numeric(names(clusx)[i]))])
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
  grp1_red <- reduction_coords[compare_list[[1]],]#; grp1_red <- rbind(grp1_red, pad_mat)
  grp2_red <- reduction_coords[compare_list[[2]],]#; grp2_red <- rbind(grp2_red, pad_mat)
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
  background_data <- reduction_coords[sample(1:nrow(reduction_coords),size=100000,replace=FALSE),]

  plt_dens_diff <- ggplot(diff12.m, aes(x = Var1, y = Var2, z=z, fill=z)) +
    geom_tile() +
    stat_contour(aes(colour = after_stat(!!str2lang("level"))), binwidth = 0.001) +
    scale_fill_gradient2(low = color_list[[1]],mid = "white",
                         high = color_list[[2]], midpoint = 0, breaks = c(minlim,maxlim),
                         labels=c(names(compare_list[1]),names(compare_list[2])),limits=c(minlim,maxlim)) +
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

  plt_leg <- ggpubr::as_ggplot(ggpubr::get_legend(p = plt_dens_diff))

  plt_dens_diff <- plt_dens_diff + theme(legend.position = "none")

  if(length(dbscan_eps_ratio)!=1) {
    dbscan_eps_ratio <- "auto"
  }
  if(length(dbscan_minPts_ratio)!=1) {
    dbscan_minPts_ratio <- "auto"
  }
  if(dbscan_eps_ratio == "auto") {
    dbscan_eps_arg <- mean(diff(dim1_range),diff(dim2_range))/20
  } else {
    dbscan_eps_arg <- mean(diff(dim1_range),diff(dim2_range))/dbscan_eps_ratio
  }
  if(dbscan_minPts_ratio == "auto") {
    dbscan_minpts_arg <- floor(nrow(background_data)/100)
  } else {
    dbscan_minpts_arg <- floor(nrow(background_data)/dbscan_minPts_ratio)
  }
  scan_islands <- dbscan::dbscan(x = background_data[,c(1:2)], eps = dbscan_eps_arg, minPts = dbscan_minpts_arg)
  background_data$island <- as.character(scan_islands$cluster)
  background_data$island_1 <- scan_islands$cluster
  override_island_col <- rep(alpha(c("black"), alpha = 1), length(unique(background_data$island)))
  names(override_island_col) <- unique(background_data$island)

  plt_dens_back <- ggplot(background_data, aes(x = dimx, y = dimy, color = island_1)) +
    geom_point(alpha = 0) +
    geom_mark_hull(concavity = 2, expand = unit(1, "mm"), aes(fill = island, filter = island != '0'),
                   color = alpha(colour = "white", alpha = 0)) +
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

  timestamp <- strftime(Sys.time(),"%Y-%m-%d_%H%M%S")
  if(add_timestamp) {
    fname_top <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                        names(compare_list)[2],"_reduction_density_difference_over_",
                        timestamp)
    fname_bottom <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_reduction_density_difference_under_",
                           timestamp)
    fname_legend <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_reduction_density_difference_legend_",
                           timestamp)
  } else {
    fname_top <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                        names(compare_list)[2],"_",tolower(reduction),
                        "_density_difference_over")
    fname_bottom <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_",tolower(reduction),
                           "_density_difference_under")
    fname_legend <- paste0(outdir,"/",names(compare_list)[1],"_vs_",
                           names(compare_list)[2],"_reduction_density_difference_legend")
  }

  ggsave(filename = paste0(fname_top,".pdf"),
         plot = plt_leg, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 900, bg = "white")
  ggsave(filename = paste0(fname_bottom,".pdf"),
         plot = plt_dens_back, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 900, bg = "white")
  ggsave(filename = paste0(fname_legend,".pdf"),
         plot = plt_dens_diff, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 900, bg = "white")
}
