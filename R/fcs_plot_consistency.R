#' @title Visualize Local Neighborhood Consistency on a 2D Embedding
#'
#' @description
#'   Plots a 2D reduction embedding with each cell coloured by a local
#'   neighbourhood consistency metric computed by
#'   \code{FCSimple::fcs_neighborhood_consistency()}.  Optionally overlays
#'   edges that are \emph{broken} (high-dimensional neighbours pulled apart)
#'   or \emph{false} (adjacencies invented by the low-dimensional embedding),
#'   directly highlighting where UMAP/t-SNE fails to preserve local structure.
#'
#' @param fcs_join_obj
#'   A list produced by \code{FCSimple::fcs_reduce_dimensions()},
#'   \code{FCSimple::fcs_cluster()}, \code{FCSimple::fcs_reduction_borders()},
#'   and \code{FCSimple::fcs_neighborhood_consistency()}.
#'
#' @param reduction
#'   Character; the 2D reduction to plot (e.g. \code{"umap"} or \code{"tsne"}).
#'   Case-insensitive.  The corresponding consistency entry (e.g.
#'   \code{"umap_2d"}) is looked up automatically.
#'
#' @param metric
#'   Character; which consistency score to colour cells by.  One of
#'   \code{"jaccard"}, \code{"continuity"}, or \code{"trustworthiness"}.
#'   Default \code{"trustworthiness"}.
#'
#' @param point_alpha
#'   Numeric; point transparency.  Default \code{0.3}.
#'
#' @param point_size
#'   Numeric; point size.  Default \code{1}.
#'
#' @param outdir
#'   Character; directory to save the plot when \code{return_plot = FALSE}.
#'   Default \code{getwd()}.
#'
#' @param figure_width, figure_height
#'   Numeric; plot dimensions in inches.  Default \code{10}.
#'
#' @param plotting_device
#'   Character; \code{"pdf"} or \code{"png"}.  Default \code{"pdf"}.
#'
#' @param add_timestamp
#'   Logical; append a timestamp to the filename.  Default \code{TRUE}.
#'
#' @param return_plot
#'   Logical; return the ggplot object instead of writing to file.
#'   Default \code{FALSE}.
#'
#' @param draw_broken_edges
#'   Logical; overlay high-dimensional neighbour pairs that were pulled apart
#'   in this embedding (coloured by \code{edge_color}).  Default \code{FALSE}.
#'
#' @param draw_false_edges
#'   Logical; overlay low-dimensional neighbour pairs that were not neighbours
#'   in high dimensions (coloured by \code{false_edge_color}).
#'   Default \code{FALSE}.
#'
#' @param edge_color
#'   Character; colour for broken (torn-apart) edges.  Default \code{"red"}.
#'
#' @param false_edge_color
#'   Character; colour for false (invented) edges.  Default \code{"blue"}.
#'
#' @param edge_alpha, edge_linewidth
#'   Numeric; styling for overlaid edges.  Defaults \code{0.5} and \code{0.3}.
#'
#' @param max_edges
#'   Integer; when drawing edges, randomly cap the number of segments to this
#'   many (per edge class) to keep the plot responsive.  Default \code{5000}.
#'
#' @param edge_seed
#'   Integer; seed for edge downsampling.  Default \code{123}.
#'
#' @param quantile_clip
#'   Numeric in (0, 1] or \code{NULL}.  When provided, the colour scale upper
#'   bound is set to this quantile of the metric so a few extreme outlier cells
#'   do not compress the informative mid-range.  Default \code{NULL} (no clip).
#'
#' @param limits
#'   Numeric vector of length 2, or \code{NULL}.  Explicit colour-scale limits.
#'   When \code{NULL} (default), the scale auto-ranges (or is quantile-clipped
#'   when \code{quantile_clip} is set).  Provide \code{c(0, 1)} to force the
#'   full theoretical range.
#'
#' @details
#'   Cell colours come from \code{fcs_join_obj$neighborhood_consistency} and
#'   are drawn with a continuous viridis scale.  Edge sets are derived from the
#'   stored high-dimensional and low-dimensional neighbour graphs in
#'   \code{fcs_join_obj$neighbors}; broken edges are those present in the
#'   high-dimensional graph but missing from the plotted embedding, while
#'   false edges are present in the embedding but absent in high dimensions.
#'
#'   \strong{Cross-embedding comparison:} the three consistency metrics are
#'   normalized to \eqn{[0, 1]} against the same high-dimensional reference and
#'   are therefore directly comparable between UMAP and t-SNE as numeric
#'   quantities.  For side-by-side visual comparison, however, pass a fixed
#'   \code{limits = c(0, 1)} (or a shared empirical range computed across the
#'   embeddings) — otherwise each plot auto-ranges independently and the same
#'   colour encodes a different value in each panel.
#'
#' @return
#'   If \code{return_plot = TRUE}, a ggplot object; otherwise the plot is
#'   written to \code{outdir} and \code{NULL} is returned invisibly.
#'
#' @examples
#' \dontrun{
#'   obj <- FCSimple::fcs_neighborhood_consistency(obj)
#'
#'   # Trustworthiness heatmap of the 2D UMAP
#'   FCSimple::fcs_plot_consistency(obj, "umap", metric = "trustworthiness")
#'
#'   # Overlay broken (red) and false (blue) edges
#'   FCSimple::fcs_plot_consistency(
#'     obj, "umap", metric = "jaccard",
#'     draw_broken_edges = TRUE, draw_false_edges = TRUE
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_neighborhood_consistency, FCSimple::fcs_reduction_borders
#'
#' @importFrom ggplot2 ggplot aes scale_color_viridis_c labs theme_bw theme
#'   ggsave geom_segment scale_color_manual
#' @importFrom ggrastr geom_point_rast
#' @importFrom viridis scale_color_viridis
#' @export
fcs_plot_consistency <- function(
    fcs_join_obj,
    reduction,
    metric = "trustworthiness",
    point_alpha = 0.3,
    point_size = 1,
    outdir = getwd(),
    figure_width = 10,
    figure_height = 10,
    plotting_device = "pdf",
    add_timestamp = TRUE,
    return_plot = FALSE,
    draw_broken_edges = FALSE,
    draw_false_edges = FALSE,
    edge_color = "red",
    false_edge_color = "blue",
    edge_alpha = 0.5,
    edge_linewidth = 0.3,
    max_edges = 5000L,
    edge_seed = 123L,
    quantile_clip = NULL,
    limits = NULL)
{
  if (!require(ggplot2, quietly = TRUE)) stop("Package 'ggplot2' is required but could not be loaded.")
  if (!require(ggrastr, quietly = TRUE)) stop("Package 'ggrastr' is required but could not be loaded.")
  if (!require(viridis, quietly = TRUE)) stop("Package 'viridis' is required but could not be loaded.")

  metric <- match.arg(metric, c("jaccard", "continuity", "trustworthiness"))

  # -- Resolve the consistency entry -----------------------------------------
  # Consistency entries are keyed by the canonical dimensional slot name
  # (e.g. "umap_2d"), mirroring $neighbors, which may differ from the
  # top-level object's legacy reduction names (e.g. "umap").
  cons <- fcs_join_obj$neighborhood_consistency
  if (is.null(cons))
    stop("No '$neighborhood_consistency' found. Run ",
         "fcs_neighborhood_consistency() first.")

  red   <- tolower(reduction)
  slot  <- paste0(red, "_2d")
  if (!slot %in% names(cons)) {
    cand <- grep(paste0("^", red, "_2d$"), names(cons), value = TRUE)
    if (length(cand) == 0L)
      stop("No neighborhood_consistency entry for '", red,
           "' (2D). Available: ", paste(names(cons), collapse = ", "))
    slot <- cand[1L]
  }

  scores <- cons[[slot]][[metric]]
  if (is.null(scores))
    stop("Metric '", metric, "' not found in neighborhood_consistency$", slot, ".")

  # -- Coordinates: canonical slot first, then legacy flat name --------------
  if (slot %in% names(fcs_join_obj)) {
    coords <- fcs_join_obj[[slot]]$coordinates
  } else if (red %in% names(fcs_join_obj)) {
    coords <- fcs_join_obj[[red]]$coordinates
  } else {
    stop("No coordinates for reduction '", red, "' found in fcs_join_obj.")
  }
  if (!is.matrix(coords))
    coords <- as.matrix(coords)

  n_cells <- nrow(coords)

  df <- data.frame(
    x     = coords[, 1L],
    y     = coords[, 2L],
    score = scores,
    stringsAsFactors = FALSE
  )
  coord_names <- colnames(coords)[1:2]

  # Determine colour-scale limits: explicit > quantile-clip > auto-range.
  if (is.null(limits)) {
    if (!is.null(quantile_clip)) {
      if (!is.numeric(quantile_clip) || length(quantile_clip) != 1L ||
          quantile_clip <= 0 || quantile_clip > 1)
        stop("'quantile_clip' must be in (0, 1].")
      lo <- min(0, min(scores, na.rm = TRUE))
      hi <- as.numeric(stats::quantile(scores, probs = quantile_clip, na.rm = TRUE))
      limits <- c(lo, hi)
    } else {
      limits <- NULL
    }
  }

  p <- ggplot(df, aes(x = x, y = y, color = score)) +
    ggrastr::geom_point_rast(alpha = point_alpha, size = point_size, stroke = 0) +
    scale_color_viridis_c(option = "inferno", limits = limits) +
    labs(x = coord_names[1L], y = coord_names[2L], color = metric)

  # -- Optional broken / false edge overlays ---------------------------------
  neighbors <- fcs_join_obj$neighbors
  if ((draw_broken_edges || draw_false_edges) &&
      !is.null(neighbors) && !is.null(neighbors$high_dim)) {
    hd <- neighbors$high_dim

    # High-dimensional reference edge list
    if (!is.null(hd$nn_idx)) {
      kh <- ncol(hd$nn_idx)
      hi <- rep(seq_len(n_cells), each = kh); hj <- as.vector(t(hd$nn_idx))
    } else {
      hi <- hd$nn_i; hj <- hd$nn_j
    }

    # Low-dimensional edge list (build from nn_idx)
    low_nb <- neighbors[[slot]]
    if (!is.null(low_nb) && !is.null(low_nb$nn_idx) && nrow(low_nb$nn_idx) == n_cells) {
      kl <- ncol(low_nb$nn_idx)
      li <- rep(seq_len(n_cells), each = kl); lj <- as.vector(t(low_nb$nn_idx))
    } else {
      li <- integer(0L); lj <- integer(0L)
    }

    if (length(hi) > 0L && length(li) > 0L) {

      # Undirected pair keys
      high_key <- paste(pmin(hi, hj), pmax(hi, hj), sep = "_")
      low_key  <- paste(pmin(li, lj), pmax(li, lj), sep = "_")

      overlay_edges <- function(from_idx, to_idx, color, n_max, s) {
        keep_undirected <- from_idx < to_idx
        fi <- from_idx[keep_undirected]; ti <- to_idx[keep_undirected]
        if (length(fi) > n_max) {
          set.seed(s)
          sel <- sample.int(length(fi), n_max)
          fi <- fi[sel]; ti <- ti[sel]
        }
        ggplot2::geom_segment(
          data = data.frame(
            x = coords[fi, 1L], y = coords[fi, 2L],
            xend = coords[ti, 1L], yend = coords[ti, 2L]),
          aes(x = x, y = y, xend = xend, yend = yend),
          color = color, alpha = edge_alpha, linewidth = edge_linewidth,
          inherit.aes = FALSE)
      }

      if (draw_broken_edges) {
        broken <- !(high_key %in% low_key)
        p <- p + overlay_edges(hi[broken], hj[broken], edge_color,
                               max_edges, edge_seed)
      }
      if (draw_false_edges) {
        false <- !(low_key %in% high_key)
        p <- p + overlay_edges(li[false], lj[false], false_edge_color,
                               max_edges, edge_seed + 1L)
      }
    } else {
      message("Neighbour-index matrices missing; skipping edge overlay.")
    }
  }

  p <- p + theme_bw(base_size = 14) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5))

  if (return_plot) {
    return(p)
  }

  if (add_timestamp) {
    fname <- paste0(outdir, "/", slot, "_", metric, "_consistency_",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"))
  } else {
    fname <- paste0(outdir, "/", slot, "_", metric, "_consistency")
  }
  if (plotting_device == "pdf") {
    ggsave(paste0(fname, ".pdf"), p, device = "pdf",
           width = figure_width, height = figure_height, units = "in", dpi = 900)
  } else if (plotting_device == "png") {
    ggsave(paste0(fname, ".png"), p, device = "png",
           width = figure_width, height = figure_height, units = "in", dpi = 900)
  } else {
    stop("plotting_device must be 'pdf' or 'png'")
  }
  message("Plot saved to: ", fname, ".", plotting_device)
  invisible(NULL)
}