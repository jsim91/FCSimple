#' @title Plot Clustered Reduction Embedding
#'
#' @description
#'   Visualizes cell clusters on a 2D reduction embedding (UMAP or tSNE) as a
#'   scatter plot with optional cluster labels, split panels, and flexible
#'   output options. Supports coloring by cluster, randomizing palettes, and
#'   saving to file or returning the ggplot object. Cluster labels can be drawn
#'   using either shadowed text or repelled text annotations.
#'
#' @param fcs_join_obj
#'   A list containing reduction coordinates and clustering results, as
#'   produced by FCSimple::fcs_reduce_dimensions() and
#'   FCSimple::fcs_cluster(). Must include:
#'   \itemize{
#'     \item `fcs_join_obj[[ tolower(reduction) ]][["coordinates"]]`: a numeric
#'       matrix of cell × 2 embedding coordinates
#'     \item `fcs_join_obj[[ tolower(algorithm) ]][["clusters"]]`: a vector of
#'       cluster IDs for each cell
#'   }
#'
#' @param algorithm
#'   Character; clustering result to visualize (e.g. `"leiden"`, `"flowsom"`).
#'
#' @param reduction
#'   Character; dimensionality‐reduction to plot. Either `"UMAP"` or `"tSNE"`.
#'
#' @param point_alpha
#'   Numeric; point transparency (alpha) for scatter (default 0.1).
#'
#' @param point_size
#'   Numeric; point size for scatter (default 1).
#'
#' @param outdir
#'   Character; directory path to save the plot when `return_plot = FALSE`
#'   (default: `getwd()`).
#'
#' @param split_factor
#'   Optional vector or factor the same length as rows of the reduction
#'   coordinates. If not `NA`, splits data by its levels and arranges panels
#'   using ggpubr::ggarrange() (default: `NA`).
#'
#' @param internal_call
#'   Logical; if `TRUE`, uses internal logic to highlight and annotate a subset
#'   of events (`keep_indices`) and clusters (`anno_indices`) rather than full scatter (default `FALSE`).
#'
#' @param anno_indices
#'   Integer vector; cluster indices to annotate when `internal_call = TRUE`
#'   (default `NULL`).
#'
#' @param keep_indices
#'   Integer vector; cell indices to highlight (in red) when
#'   `internal_call = TRUE` (default `NA`).
#'
#' @param figure_width
#'   Numeric; width in inches for saving plot (default 10).
#'
#' @param figure_height
#'   Numeric; height in inches for saving plot (default 10).
#'
#' @param plotting_device
#'   Character; `"pdf"` or `"png"` to select output device when writing file
#'   (default `"pdf"`).
#'
#' @param annotate_text_size
#'   Numeric; font size for cluster labels (use `NA` to disable, default 5).
#'
#' @param annotation_method
#'   Character; method for drawing cluster labels. Options are:
#'   \itemize{
#'     \item `"shadowtext"` (default): labels with shadowed outlines using
#'       **shadowtext::geom_shadowtext**
#'     \item `"repel"`: labels with repelling force to avoid overlap using
#'       **ggrepel::geom_text_repel**
#'   }
#'
#' @param title_size
#'   Numeric; font size for plot title (default 14).
#'
#' @param return_plot
#'   Logical; if `TRUE` (default), returns the ggplot2 object; if `FALSE`,
#'   writes the plot to file and invisibly returns `NULL`.
#'
#' @param randomize_colors
#'   Logical; if `TRUE`, shuffles cluster‐color assignment (default `FALSE`).
#'
#' @param color_random_seed
#'   Integer; seed for random color assignment when `randomize_colors = TRUE`
#'   (default 123).
#'
#' @param color_clusters
#'   Logical; if `TRUE` (default), color points by cluster; if `FALSE`, all
#'   points are drawn in a single color.
#'
#' @param force_title
#'   Logical; if `TRUE`, forces display of the reduction name as a title even
#'   in split panels (default `FALSE`).
#'
#' @param sample_equally
#'   Logical; if `TRUE`, downsamples each split‐group to equal size before plotting (default `TRUE`).
#'
#' @param cluster_substitute_names
#'   Named character vector; optional mapping from original cluster IDs to
#'   replacement labels. The names of this vector must exactly match the unique
#'   cluster IDs in `fcs_join_obj[[ tolower(algorithm) ]][["clusters"]]`;
#'   a mismatch will trigger an error. If supplied (i.e. not `NA`), cluster
#'   numbers in the plot are replaced by their mapped labels. Labels may
#'   include the newline characters for multi‐line annotations.
#'   Default is `NA` (no substitution).
#'
#' @param add_timestamp
#'   Logical; if `TRUE` (default), appends a timestamp (`_YYYY-MM-DD_HHMMSS`)
#'   to filenames when saving.
#'
#' @details
#' \enumerate{
#'   \item Extracts embedding coordinates and cluster IDs.
#'   \item Builds a data.frame for ggplot, mapping clusters to colors via a
#'      default HCL palette (or randomized if requested).
#'   \item Computes median centroids for each cluster for annotation.
#'   \item If `split_factor` is `NA`, draws a single scatter plot. Otherwise splits
#'      cells by factor levels, optionally downsamples equally, and arranges
#'      subplots with ggpubr::ggarrange().
#'   \item If `internal_call = TRUE`, highlights cells in `keep_indices` and
#'      annotates clusters in `anno_indices`.
#'   \item Cluster labels are drawn using the method specified in
#'      `annotation_method` (`"shadowtext"` or `"repel"`).
#'   \item When `return_plot = FALSE`, saves to `outdir` as
#'      `<algorithm>_<reduction>_labeled[(_timestamp)].pdf` or `.png`.
#' }
#'
#' @return
#'   If `return_plot = TRUE`, a ggplot2 object (or ggarrange object) is returned.
#'   If `return_plot = FALSE`, the plot is written to file and the function
#'   invisibly returns `NULL`.
#'
#' @examples
#' \dontrun{
#'   files   <- list.files("data/", "\\.fcs$", full.names = TRUE)
#'   joined  <- FCSimple::fcs_join(files)
#'   reduced <- FCSimple::fcs_reduce_dimensions(joined, method = "UMAP")
#'   clustered <- FCSimple::fcs_cluster(reduced, algorithm = "leiden")
#'
#'   # Return the plot object
#'   p <- FCSimple::fcs_plot_reduction(clustered, "leiden", "UMAP")
#'   print(p)
#'
#'   # Save without returning, as PNG, no timestamp
#'   FCSimple::fcs_plot_reduction(
#'     clustered, "leiden", "UMAP",
#'     return_plot      = FALSE,
#'     plotting_device  = "png",
#'     figure_width     = 8,
#'     figure_height    = 8,
#'     add_timestamp    = FALSE
#'   )
#'
#'   # Split panels by sample
#'   p2 <- FCSimple::fcs_plot_reduction(
#'     clustered, "leiden", "UMAP",
#'     split_factor  = clustered$source,
#'     sample_equally = FALSE
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_reduce_dimensions, FCSimple::fcs_cluster,
#'   ggplot2::ggplot, ggrastr::geom_point_rast, ggpubr::ggarrange,
#'   shadowtext::geom_shadowtext, ggrepel::geom_text_repel
#'
#' @importFrom grDevices hcl
#' @importFrom ggplot2 ggplot aes scale_color_manual labs theme_void theme ggtitle annotate ggsave
#' @importFrom ggrastr geom_point_rast
#' @importFrom ggpubr ggarrange
#' @importFrom shadowtext geom_shadowtext
#' @export
fcs_plot_reduction <- function(fcs_join_obj, algorithm, reduction, point_alpha = 0.1, point_size = 1, outdir = getwd(),
                               split_factor = NA, internal_call = FALSE, anno_indices = NULL, keep_indices = NA,
                               figure_width = 10, figure_height = 10, plotting_device = "pdf", annotate_text_size = 5,
                               title_size = 14, return_plot = TRUE, randomize_colors = FALSE, color_random_seed = 123,
                               color_clusters = TRUE, force_title = FALSE, sample_equally = TRUE,
                               cluster_substitute_names = NA, add_timestamp = TRUE, annotation_method = 'shadowtext')
{
  # use annotate_text_size = NA to produce a plot without cluster annotations
  require(ggplot2)
  require(ggrastr)
  require(ggpubr)
  require(shadowtext)

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  cluster_numbers <- as.character(fcs_join_obj[[tolower(algorithm)]][["clusters"]])
  if(!is.na(cluster_substitute_names[1])) {
    if(mean(names(cluster_substitute_names) %in% unique(cluster_numbers))!=1) {
      stop("error in argument 'cluster_substitute_names': length of 'cluster_substitute_names' does not match number of clusters.")
    } else {
      print("Substituted cluster names will replace cluster numbers. You may use '\n' character to create multi-line annotations. 'My\nCluster' will put 'My' on first line and 'Cluster' on second line.")
      cluster_numbers <- as.character(cluster_substitute_names[cluster_numbers])
    }
  }
  uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]

  # https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  colorby <- gg_color_hue(n = length(uclus))
  if(randomize_colors) {
    set.seed(color_random_seed)
    names(colorby) <- sample(x = uclus, size = length(uclus), replace = FALSE)
  } else {
    names(colorby) <- uclus
  }

  xval <- rep(NA,times=length(uclus)); names(xval) <- uclus; yval <- xval
  for(i in 1:length(xval)) {
    xval[i] <- median(reduction_coords[,1][which(cluster_numbers==names(xval)[i])])
    yval[i] <- median(reduction_coords[,2][which(cluster_numbers==names(yval)[i])])
  }

  plt_input <- cbind(reduction_coords,data.frame(cluster = cluster_numbers))
  plt_input$cluster <- factor(plt_input$cluster)
  if(!internal_call) {
    pl_fun <- function(plin, ptalpha = point_alpha, xanno = xval,
                       yanno = yval, sizeanno = annotate_text_size, ftitle = force_title,
                       color_clus = color_clusters, ptsize = point_size)
    {
      cnamex <- colnames(plin)[1]; cnamey <- colnames(plin)[2]
      colnames(plin)[1:2] <- c("valx","valy") # bypass aes_string, which is now deprecated
      if(color_clus) {
        mypl <- ggplot(data = plin, mapping = aes(x = valx,
                                                  y = valy,
                                                  color = cluster)) +
          ggrastr::geom_point_rast(alpha = ptalpha, size = ptsize) +
          scale_color_manual(values = colorby) +
          labs(x = cnamex, y = cnamey)
      } else {
        mypl <- ggplot(data = plin, mapping = aes(x = valx,
                                                  y = valy)) +
          ggrastr::geom_point_rast(alpha = ptalpha, size = ptsize) +
          labs(x = cnamex, y = cnamey)
      }
      if(!is.na(sizeanno)) {
        mypl <- mypl +
          annotate("shadowtext", x = xanno, y = yanno, label = names(xanno), size = sizeanno)
      }
      mypl <- mypl + theme_void() +
        theme(legend.position = "none")
      if(ftitle) {
        mypl <- mypl + ggtitle(colnames(plin)[ncol(plin)]) + theme(plot.title = element_text(hjust = 0.5, size = title_size))
      }
      return(mypl)
    }
    if(is.na(split_factor[1])) {
      plt_reduction <- pl_fun(plin = plt_input)
    } else {
      # split the plot by split_factor where split_factor is a vector of length == nrow(reduction) that gives identity to the reduction rows
      split_reduction <- split(x = plt_input, f = factor(split_factor))
      for(i in 1:length(split_reduction)) {
        split_reduction[[i]] <- cbind(as.data.frame(split_reduction[[i]]), data.frame(var1 = rep(1,nrow(split_reduction[[i]]))))
        colnames(split_reduction[[i]])[ncol(split_reduction[[i]])] <- names(split_reduction)[i]
      }
      if(sample_equally) {
        sample_size <- min(as.numeric(sapply(split_reduction,nrow)))
        for(i in 1:length(split_reduction)) {
          if(nrow(split_reduction[[i]])>sample_size) {
            set.seed(123)
            split_reduction[[i]] <- split_reduction[[i]][sample(1:nrow(split_reduction[[i]]),sample_size,replace=F),]
          }
        }
      }
      outplots <- lapply(X = split_reduction, pl_fun, ftitle = force_title)
      plotnrow <- ifelse(length(outplots)<=4,1,floor(sqrt(length(outplots))))
      plotncol <- ifelse(length(outplots)<=4,length(outplots),ceiling(sqrt(length(outplots))))
      plt_reduction <- ggpubr::ggarrange(plotlist = outplots, nrow = plotnrow, ncol = plotncol)
    }
    if(add_timestamp) {
      fname <- paste0(outdir,"/",tolower(algorithm),"_",tolower(reduction),"_labeled_",
                      strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
    } else {
      fname <- paste0(outdir,"/",tolower(algorithm),"_",tolower(reduction),"_labeled")
    }
  } else {
    plt_reduction <- ggplot(data = plt_input[-keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                       y = colnames(reduction_coords)[2])) +
      ggrastr::geom_point_rast(alpha = point_alpha, color = "grey") +
      ggrastr::geom_point_rast(data = plt_input[keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                        y = colnames(reduction_coords)[2]),
                               alpha = point_alpha, color = "red")
    if(annotation_method=='shadowtext') {
      plt_reduction <- plt_reduction + annotate("shadowtext", x = xval, y = yval, label = names(xval), size = annotate_text_size)
    } else if(annotation_method=='repel') {
      require(ggrepel)
      color_text_add <- data.frame(UMAP1 = xval, UMAP2 = yval, cluster = names(xval))
      plt_reduction <- plt_reduction + ggrepel::geom_text_repel(data = color_text_add,
                                                                mapping = aes(x = UMAP1, y = UMAP2, label = cluster), color = "white",
                                                                size = annotate_text_size, bg.color = "black", bg.r = 0.04, seed = 123)
    }
      plt_reduction <- theme_void() +
      theme(legend.position = "none")
    if(add_timestamp) {
      fname <- paste0(outdir,"/islands_selected_for_by_dbscan_",
                      strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
    } else {
      fname <- paste0(outdir,"/islands_selected_for_by_dbscan_")
    }
  }
  if(return_plot) {
    return(plt_reduction)
  } else {
    if(plotting_device=="pdf") {
      ggsave(filename = paste0(fname,".pdf"),
             plot = plt_reduction, device = "pdf", width = figure_width, height = figure_height,
             units = "in", dpi = 900, bg = "white")
    } else if(plotting_device=="png"){
      ggsave(filename = paste0(fname,".png"),
             plot = plt_reduction, device = "png", width = figure_width, height = figure_height,
             units = "in", dpi = 900, bg = "white")
    }
  }
}
