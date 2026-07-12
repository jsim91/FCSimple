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
#'   Logical; if `TRUE`, returns the ggplot2 object; if `FALSE` (default),
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
#'   reduced <- FCSimple::fcs_reduce_dimensions(joined, algorithm = "umap")
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
                               figure_width = 10, figure_height = 10, plotting_device = "pdf", annotate_text_size = 8,
                               title_size = 14, return_plot = FALSE, randomize_colors = FALSE, color_random_seed = 123,
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

  # Helper: generate a perceptually-balanced categorical palette matching
  # sc0rch's approach — RColorBrewer Set2 (≤8), Paired (≤12), then
  # farthest-point sampling in CIE Lab for larger n.
  fcs_color_palette <- function(n) {
    if (n <= 8L)  return(RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)])
    if (n <= 12L) return(RColorBrewer::brewer.pal(max(3L, n), "Paired")[seq_len(n)])

    # Farthest-point sampling in CIE Lab over an HCL candidate grid
    h_vals <- seq(0, 355, by = 5)
    c_vals <- seq(40, 100, by = 10)
    l_vals <- seq(40, 85,  by = 5)
    grid   <- expand.grid(H = h_vals, C = c_vals, L = l_vals)
    hex    <- grDevices::hcl(h = grid$H, c = grid$C, l = grid$L, fixup = FALSE)
    ok     <- !is.na(hex)
    hex    <- hex[ok]
    if (length(hex) < n) return(scales::hue_pal()(n))

    rgb_m <- t(grDevices::col2rgb(hex)) / 255
    lab   <- grDevices::convertColor(rgb_m, from = "sRGB", to = "Lab")

    sel      <- integer(n)
    sel[1L]  <- which.max(lab[, 2L]^2 + lab[, 3L]^2)   # seed: most chromatic
    min_d    <- sqrt(rowSums(sweep(lab, 2, lab[sel[1L], ])^2))
    for (i in seq.int(2L, n)) {
      sel[i] <- which.max(min_d)
      min_d  <- pmin(min_d, sqrt(rowSums(sweep(lab, 2, lab[sel[i], ])^2)))
    }
    hex[sel]
  }

  # Helper: assign colors so spatially nearby clusters get maximally different
  # colors via a DSatur (degree-of-saturation) greedy algorithm using full
  # CIE Lab perceptual distance.
  #   1) Build a k-nearest-neighbor graph (k=4) on cluster centroids.
  #   2) At each step pick the uncolored cluster with the most already-colored
  #      neighbors (DSatur) — breaking ties by the minimum spatial distance to
  #      any already-colored cluster.
  #   3) Assign the unused palette color whose minimum CIE Lab perceptual
  #      distance to those already-colored neighbors is greatest.
  optimize_cluster_colors <- function(xval, yval, base_colors) {
    n_clusters   <- length(xval)
    n_colors     <- length(base_colors)
    cluster_names <- names(xval)

    # Euclidean distance matrix between cluster centroids
    dist_matrix <- as.matrix(stats::dist(cbind(xval, yval), method = "euclidean"))

    # Build k-NN graph on centroids (k = 4, or fewer for small n)
    k <- min(4L, n_clusters - 1L)
    neighbors <- t(apply(dist_matrix, 1L, function(d) order(d)[seq_len(k + 1L)][-1L]))
    rownames(neighbors) <- cluster_names

    # Convert base_colors to CIE Lab once for perceptual distance calculations
    rgb_m    <- t(grDevices::col2rgb(base_colors)) / 255
    lab_cols <- grDevices::convertColor(rgb_m, from = "sRGB", to = "Lab")
    rownames(lab_cols) <- seq_len(n_colors)

    color_assignment <- integer(n_clusters)
    names(color_assignment) <- cluster_names
    color_used <- logical(n_colors)

    # Seed: cluster closest to the global centroid
    start_cluster <- names(which.min(sqrt((xval - mean(xval))^2 + (yval - mean(yval))^2)))
    color_assignment[start_cluster] <- 1L
    color_used[1L] <- TRUE
    colored_set <- start_cluster

    while (length(colored_set) < n_clusters) {
      uncolored <- setdiff(cluster_names, colored_set)

      # DSatur: count already-colored neighbors for each uncolored cluster
      saturation <- sapply(uncolored, function(uc) {
        sum(neighbors[uc, ] %in% colored_set)
      })

      # Break ties by minimum distance to any already-colored cluster
      min_dist_to_colored <- sapply(uncolored, function(uc) {
        min(dist_matrix[uc, colored_set])
      })

      # Pick the uncolored cluster with most colored neighbors;
      # ties broken by smallest minimum distance to a colored cluster
      best_idx <- which.max(saturation - min_dist_to_colored / max(dist_matrix))
      next_cluster <- uncolored[best_idx]

      # Gather colors already used by this cluster's already-colored neighbors
      nbr_colored <- intersect(neighbors[next_cluster, ], colored_set)
      nbr_col_idx <- color_assignment[nbr_colored]

      available <- which(!color_used)

      if (length(available) > 0L && length(nbr_col_idx) > 0L) {
        # Choose the available color whose minimum Lab distance
        # to any neighbor color is maximized
        scores <- sapply(available, function(ci) {
          min(sapply(nbr_col_idx, function(ni) {
            sqrt(sum((lab_cols[ci, ] - lab_cols[ni, ])^2))
          }))
        })
        best_color <- available[which.max(scores)]
      } else if (length(available) > 0L) {
        best_color <- available[1L]
      } else {
        # All colors used — pick the first color (edge case)
        best_color <- 1L
      }

      color_assignment[next_cluster] <- best_color
      color_used[best_color] <- TRUE
      colored_set <- c(colored_set, next_cluster)
    }

    result        <- base_colors[color_assignment]
    names(result) <- cluster_names
    result
  }

  # Calculate cluster centers (median for color optimization, mean for annotation)
  xval <- rep(NA, length(uclus)); names(xval) <- uclus
  yval <- xval
  xmean <- xval; ymean <- xval
  for (i in seq_along(xval)) {
    xval[i]  <- median(reduction_coords[, 1][cluster_numbers == names(xval)[i]])
    yval[i]  <- median(reduction_coords[, 2][cluster_numbers == names(yval)[i]])
    xmean[i] <- mean(reduction_coords[, 1][cluster_numbers == names(xmean)[i]])
    ymean[i] <- mean(reduction_coords[, 2][cluster_numbers == names(ymean)[i]])
  }

  # Assign cluster colors
  base_colors <- fcs_color_palette(length(uclus))
  if (randomize_colors) {
    set.seed(color_random_seed)
    colorby        <- base_colors
    names(colorby) <- sample(uclus, length(uclus), replace = FALSE)
  } else {
    colorby <- optimize_cluster_colors(xval, yval, base_colors)
  }

  plt_input <- cbind(reduction_coords,data.frame(cluster = cluster_numbers))
  plt_input$cluster <- factor(plt_input$cluster)
  if(!internal_call) {
    pl_fun <- function(plin, ptalpha = point_alpha, xanno = xmean, ameth = annotation_method,
                       yanno = ymean, sizeanno = annotate_text_size, ftitle = force_title,
                       color_clus = color_clusters, ptsize = point_size)
    {
      cnamex <- colnames(plin)[1]; cnamey <- colnames(plin)[2]
      colnames(plin)[1:2] <- c("valx","valy") # bypass aes_string, which is now deprecated
      if(color_clus) {
        mypl <- ggplot(data = plin, mapping = aes(x = valx,
                                                  y = valy,
                                                  color = cluster)) +
          ggrastr::geom_point_rast(alpha = ptalpha, size = ptsize, stroke = 0) +
          scale_color_manual(values = colorby) +
          labs(x = cnamex, y = cnamey)
      } else {
        mypl <- ggplot(data = plin, mapping = aes(x = valx,
                                                  y = valy)) +
          ggrastr::geom_point_rast(alpha = ptalpha, size = ptsize, stroke = 0) +
          labs(x = cnamex, y = cnamey)
      }
      if(!is.na(sizeanno) && !is.na(ameth)) {
        if(ameth=='shadowtext') {
          centroids <- data.frame(valx = xanno, valy = yanno,
                                  cluster = names(xanno), stringsAsFactors = FALSE)
          mypl <- mypl + shadowtext::geom_shadowtext(
            data        = centroids,
            mapping     = aes(x = valx, y = valy, label = cluster),
            color       = "white",
            bg.color    = "black",
            bg.r        = 0.15,
            size        = sizeanno,
            fontface    = "bold",
            inherit.aes = FALSE
          )
        } else if(ameth=='repel') {
          require(ggrepel)
          color_text_add <- data.frame(valx = xanno, valy = yanno, cluster = names(xanno))
          mypl <- mypl + ggrepel::geom_text_repel(data = color_text_add,
                                                  mapping = aes(x = valx, y = valy, label = cluster),
                                                  color = "white",
                                                  size = sizeanno,
                                                  fontface = "bold",
                                                  box.padding = 0.3,
                                                  point.padding = 0.2,
                                                  segment.alpha = 0.4,
                                                  inherit.aes = FALSE)
        }
      }
      mypl <- mypl + theme_bw(base_size = title_size) +
        theme(legend.position = "none",
              plot.title = element_text(face = "bold"))
      if(ftitle) {
        mypl <- mypl + ggtitle(colnames(plin)[ncol(plin)]) + 
        theme(plot.title = element_text(hjust = 0.5, size = title_size))
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
      ggrastr::geom_point_rast(alpha = point_alpha, color = "grey", stroke = 0) +
      ggrastr::geom_point_rast(data = plt_input[keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                        y = colnames(reduction_coords)[2]),
                               alpha = point_alpha, color = "red", stroke = 0)
    if(!is.na(annotation_method) && annotation_method=='shadowtext') {
      centroids <- data.frame(x = xmean, y = ymean, cluster = names(xmean), stringsAsFactors = FALSE)
      plt_reduction <- plt_reduction + shadowtext::geom_shadowtext(
        data        = centroids,
        mapping     = aes(x = x, y = y, label = cluster),
        color       = "white",
        bg.color    = "black",
        bg.r        = 0.15,
        size        = annotate_text_size,
        fontface    = "bold",
        inherit.aes = FALSE
      )
    } else if(!is.na(annotation_method) && annotation_method=='repel') {
      require(ggrepel)
      color_text_add <- data.frame(UMAP1 = xmean, UMAP2 = ymean, cluster = names(xmean))
      plt_reduction <- plt_reduction + ggrepel::geom_text_repel(data = color_text_add, force = 0, force_pull = Inf,
                                                                mapping = aes(x = UMAP1, y = UMAP2, label = cluster), color = "white",
                                                                size = annotate_text_size,
                                                                fontface = "bold",
                                                                bg.color = "black", bg.r = 0.15, seed = 123)
    }
      plt_reduction <- plt_reduction + theme_bw(base_size = 22) +
      theme(legend.position = "none",
            plot.title = element_text(face = "bold"))
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
    # ggsave cannot handle ggarrange objects; use open-device + print for those
    if(inherits(plt_reduction, "ggarrange")) {
      if(plotting_device=="pdf") {
        pdf(file = paste0(fname,".pdf"), width = figure_width, height = figure_height)
        print(plt_reduction)
        dev.off()
      } else if(plotting_device=="png") {
        png(filename = paste0(fname,".png"), width = figure_width, height = figure_height,
            units = "in", res = 900)
        print(plt_reduction)
        dev.off()
      }
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
}