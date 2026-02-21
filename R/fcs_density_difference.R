#' @title Plot Difference in Cell Density on a 2D Reduction Embedding
#'
#' @description
#'   Computes and plots the difference between two sample‐group density estimates
#'   on a 2D reduction (UMAP or tSNE). Generates filled‐contour difference maps,
#'   underlying enrichment hulls, and a combined figure, then saves each as a PDF.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join() and
#'   FCSimple::fcs_reduce_dimensions(), containing:
#'   - a named element “umap” or “tsne” with component “coordinates” (cell × 2 matrix)
#'   - clustering results under one of “leiden”, “louvain”, “flowsom”,
#'     or “phenograph”
#'   - a “source” vector of sample identifiers matching rows of “coordinates”
#'
#' @param reduction
#'   Character; which embedding to use. Either “UMAP” (default) or “tSNE”.
#'
#' @param compare_list
#'   Named list of length 2, each element a character vector of sample IDs
#'   to compare. Names must be unique and identical between `compare_list` and `color_list`.
#'
#' @param color_list
#'   Named list of length 2, each element a single color (e.g. “red”). Names must
#'   match those in `compare_list` and determine the low/high fill endpoints.
#'
#' @param n_kde
#'   Integer; number of grid points per axis for `MASS::kde2d()` (default 200).
#'
#' @param outdir
#'   Character; directory path where PDFs will be saved (default `getwd()`).
#'
#' @param axis_title_text_size
#'   Numeric; font size for axis titles (default 12).
#'
#' @param dbscan_eps
#'   Numeric or “auto”; epsilon parameter for `dbscan::dbscan()` on background
#'   points. If “auto”, set to mean axis‐ranges/20 (default “auto”).
#'
#' @param dbscan_minPts
#'   Numeric or “auto”; minPts for `dbscan::dbscan()`. If “auto”, set to
#'   floor(nrow(background)/100) (default “auto”).
#'
#' @param legend_label_text_size
#'   Numeric; font size for legend labels (default 12).
#'
#' @param annotate_clusters
#'   Logical; if `TRUE` (default), cluster IDs are annotated at each cluster’s
#'   median coordinate.
#'
#' @param cluster_algorithm
#'   Character; which clustering result to annotate. One of
#'   “leiden”, “louvain”, “flowsom”, or “phenograph” (default first).
#'
#' @param cluster_number_annotation_size
#'   Numeric; font size for cluster‐ID annotations (default 5).
#'
#' @param contour_bin_width
#'   Numeric or `NULL`; bin-width for `stat_contour()`. If `NULL` (default),
#'   computed automatically as `diff(range(density_difference)) / 20`.
#'
#' @param legend_orientation
#'   Character; “horizontal” or “vertical” (default “horizontal”) orientation
#'   for the difference–fill legend.
#'
#' @param figure_width
#'   Numeric; width (in inches) of individual PDF files (default 8).
#'
#' @param figure_height
#'   Numeric; height (in inches) of individual PDF files (default 8).
#'
#' @param hull_radius
#'   Numeric; radius (in mm) expansion around DBSCAN clusters for hull plotting
#'   (default 1.5).
#'
#' @param add_timestamp
#'   Logical; if `TRUE` (default), appends timestamp to PDF filenames.
#'
#' @param hull_concavity
#'   Numeric; concavity parameter for hulls drawn by `ggforce::geom_mark_hull()`
#'   (default 2).
#'
#' @details
#'   1. Splits cells into two groups via `compare_list` and downsamples equally.  
#'   2. Computes 2D kernel density estimates (`MASS::kde2d`) on each group.  
#'   3. Calculates difference (`group2 – group1`) and plots with `geom_tile()` +
#'      `stat_contour()`.  
#'   4. Computes background islands via `dbscan::dbscan()` and draws enrichment
#'      hulls with `geom_mark_hull()`.  
#'   5. Arranges difference plot, hull plot, and shared legend side‐by‐side.  
#'   6. Saves four PDFs to `outdir`:  
#'      - `<grp1>_vs_<grp2>_<reduction>_density_difference_under.pdf`  
#'      - `<grp1>_vs_<grp2>_<reduction>_density_difference_over.pdf`  
#'      - `<grp1>_vs_<grp2>_<reduction>_density_difference_legend.pdf`  
#'      - `<grp1>_vs_<grp2>_<reduction>_density_difference_combined.pdf`  
#'
#' @return
#'   Invisibly returns `NULL` after saving the PDF files.
#'
#' @examples
#' \dontrun{
#'   # Compute reduction first
#'   obj <- FCSimple::fcs_reduce_dimensions(joined, method = "UMAP")
#'
#'   # Define groups and colors
#'   compare <- list(Healthy = c("S1","S2"), Diseased = c("S3","S4"))
#'   colors  <- list(Healthy = "blue", Diseased = "red")
#'
#'   # Plot and save density‐difference maps
#'   FCSimple::fcs_plot_reduction_difference(
#'     obj,
#'     reduction = "UMAP",
#'     compare_list = compare,
#'     color_list   = colors,
#'     n_kde        = 150,
#'     outdir       = "~/results",
#'     annotate_clusters = FALSE
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_plot_reduction_density,
#'   FCSimple::fcs_reduce_dimensions,
#'   MASS::kde2d,
#'   dbscan::dbscan,
#'   ggforce::geom_mark_hull
#'
#' @importFrom MASS kde2d
#' @importFrom ggplot2 ggplot aes geom_tile stat_contour scale_fill_gradient2 scale_colour_gradient2 coord_cartesian guides labs annotate theme_bw theme ggsave
#' @importFrom ggforce geom_mark_hull
#' @importFrom dbscan dbscan
#' @importFrom ggpubr ggarrange get_legend
#' @importFrom circlize colorRamp2
#' @importFrom scales muted
#' @importFrom grid unit grid.grabExpr
#' @export
fcs_plot_reduction_difference <- function(fcs_join_obj, reduction = c("UMAP","tSNE"),
                                          compare_list, color_list, n_kde = 200,
                                          outdir = getwd(), axis_title_text_size = 12,
                                          dbscan_eps = "auto", dbscan_minPts = "auto",
                                          legend_label_text_size = 12, annotate_clusters = TRUE,
                                          cluster_algorithm = c("leiden","flowsom","louvain","phenograph"),
                                          cluster_number_annotation_size = 5, contour_bin_width = NULL,
                                          legend_orientation = c("horizontal","vertical"),
                                          figure_width = 8, figure_height = 8, hull_radius = 1.5,
                                          add_timestamp = TRUE, hull_concavity = 2)
{
  require(ggplot2)
  require(MASS)
  require(scales)
  require(ggrastr)
  require(ggforce); require(concaveman)
  require(dbscan)
  require(ggpubr)
  require(shadowtext)
  require(ComplexHeatmap)

  reduction        <- match.arg(reduction)
  legend_orientation <- match.arg(legend_orientation)
  cluster_algorithm  <- match.arg(cluster_algorithm)

  if (!dir.exists(outdir)) {
    stop("'outdir' does not exist: ", outdir)
  }

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  reduction_names <- colnames(reduction_coords)[1:2]
  colnames(reduction_coords)[1:2] <- c("dimx","dimy")
  dim1_range <- range(reduction_coords[,1]); dim2_range <- range(reduction_coords[,2])

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
  grp1_red <- reduction_coords[grp1_reduction_loc, ]
  grp2_red <- reduction_coords[grp2_reduction_loc, ]
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
  diff12 <- grp1_kde2d
  diff12$z <- grp2_kde2d$z - grp1_kde2d$z
  diff12.m <- data.frame(
    Var1 = rep(diff12$x, times = n_kde),
    Var2 = rep(diff12$y, each  = n_kde),
    z    = as.vector(diff12$z)
  )

  if (is.null(contour_bin_width)) {
    contour_bin_width <- diff(range(diff12.m$z)) / 20
  }

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
  plt_dens_diff <- plt_dens_diff +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction= legend_orientation,
          legend.title = element_blank(),
          legend.text = element_text(size = legend_label_text_size, hjust = 0.5, vjust = 1),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(face = "bold", size = axis_title_text_size))
  if (annotate_clusters) {
    plt_dens_diff <- plt_dens_diff +
      annotate("shadowtext", x = clusx, y = clusy, label = names(clusx), size = cluster_number_annotation_size)
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
  lgd_export <- if (legend_orientation == "horizontal") lgd_h else lgd_v
  lg_arr <- ggpubr::as_ggplot(grid.grabExpr(expr = draw(lgd_export, x = unit(0.5, "npc"), y = unit(0.5, "npc"))))

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

  plt_dens_back <- ggplot(background_data, aes(x = dimx, y = dimy)) +
    geom_mark_hull(concavity = hull_concavity, expand = unit(0.25, "mm"),
                   aes(fill = island, filter = island != '0'),
                   color = alpha(colour = "white", alpha = 0), radius = unit(hull_radius, "mm")) +
    scale_fill_manual(values = override_island_col) +
    coord_cartesian(xlim = dim1_range, ylim = dim2_range, expand = FALSE) +
    guides(fill = "none", color = "none") +
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

  ts_suffix <- if (add_timestamp) paste0("_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")) else ""
  base <- paste0(names(compare_list)[1], "_vs_", names(compare_list)[2], "_", tolower(reduction))
  fname_top    <- file.path(outdir, paste0(base, "_density_difference_over",   ts_suffix))
  fname_bottom <- file.path(outdir, paste0(base, "_density_difference_under",  ts_suffix))
  fname_legend <- file.path(outdir, paste0(base, "_density_difference_legend", ts_suffix))
  fname        <- file.path(outdir, paste0(base, "_density_difference",        ts_suffix))

  rm_box <- theme(plot.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_blank())

  umap_arr <- ggpubr::ggarrange(plotlist = list(plt_dens_back + rm_box, plt_dens_diff + rm_box + theme(legend.position = "none", axis.title = element_blank())), nrow = 1)

  ggsave(filename = paste0(fname_bottom,".pdf"),
         plot = plt_dens_back, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 300, bg = "white")
  ggsave(filename = paste0(fname_top,".pdf"),
         plot = plt_dens_diff + theme(legend.position = "none"), device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 300, bg = "white")
  ggsave(filename = paste0(fname_legend,".pdf"),
         plot = lg_arr, device = "pdf",
         width = figure_width, height = figure_height,
         units = "in", dpi = 300, bg = "white")
  ggsave(filename = paste0(fname,".pdf"),
         plot = umap_arr, device = "pdf",
         width = figure_width*2, height = figure_height,
         units = "in", dpi = 300, bg = "white")
}
