#' @title Plot Border Points on a 2D Reduction Embedding
#'
#' @description
#'   Visualizes border points — cells identified as lying at cluster boundaries
#'   or cluster peripheries — on a 2D reduction embedding (UMAP or t-SNE).
#'   Non-border points are drawn as a rasterized background layer coloured by
#'   cluster membership; border points are drawn on top as crisp vector points
#'   in black for maximum visibility.
#'
#' @param fcs_join_obj
#'   A list containing reduction coordinates, clustering results, and border
#'   scores, as produced by \code{FCSimple::fcs_reduce_dimensions()},
#'   \code{FCSimple::fcs_cluster()}, and
#'   \code{FCSimple::fcs_reduction_borders()}.
#'
#' @param algorithm
#'   Character; clustering result to colour non-border points by
#'   (e.g. \code{"leiden"}, \code{"flowsom"}).  Also used to label
#'   the plot title.  Case-insensitive.
#'
#' @param reduction
#'   Character; dimensionality-reduction to plot.  Either \code{"UMAP"} or
#'   \code{"tSNE"} (2D only).  Case-insensitive.
#'
#' @param borders_slot
#'   Character; which border-score data.frame to use from
#'   \code{fcs_join_obj$reduction_borders} (e.g. \code{"umap_2d"},
#'   \code{"umap_3d"}, \code{"tsne_2d"}, \code{"high_dim"}).
#'   Case-insensitive.
#'
#' @param border_threshold
#'   Numeric; cells with \code{border_score >= border_threshold} are
#'   drawn as border points.  Default \code{0.3}.
#'
#' @param point_alpha
#'   Numeric; transparency for non-border points (0–1).  Default \code{0.1}.
#'
#' @param point_size
#'   Numeric; base point size.  Default \code{1}.
#'
#' @param border_size_factor
#'   Numeric; multiplier for border point size relative to
#'   \code{point_size}.  Default \code{3}.
#'
#' @param outdir
#'   Character; directory to save the plot when \code{return_plot = FALSE}.
#'   Default \code{getwd()}.
#'
#' @param figure_width
#'   Numeric; plot width in inches.  Default \code{10}.
#'
#' @param figure_height
#'   Numeric; plot height in inches.  Default \code{10}.
#'
#' @param plotting_device
#'   Character; output format: \code{"pdf"} or \code{"png"}.
#'   Default \code{"pdf"}.
#'
#' @param annotate_text_size
#'   Numeric; font size for cluster centroid labels.  Set \code{NA} to
#'   hide labels.  Default \code{8}.
#'
#' @param annotation_method
#'   Character; \code{"shadowtext"} (default) or \code{"repel"} for
#'   centroid label placement.
#'
#' @param title_size
#'   Numeric; font size for the plot title.  Default \code{14}.
#'
#' @param return_plot
#'   Logical; if \code{TRUE}, returns the ggplot object.  If \code{FALSE}
#'   (default), writes to file and invisibly returns \code{NULL}.
#'
#' @param randomize_colors
#'   Logical; if \code{TRUE}, shuffles cluster-colour assignment.
#'   Default \code{FALSE}.
#'
#' @param color_random_seed
#'   Integer; seed for colour randomisation.  Default \code{123}.
#'
#' @param cluster_substitute_names
#'   Named character vector; optional cluster ID → label mapping.
#'   Default \code{NA}.
#'
#' @param add_timestamp
#'   Logical; if \code{TRUE} (default), appends a timestamp to the filename.
#'
#' @param draw_edges
#'   Logical; if \code{TRUE}, border points are connected by edges tracing
#'   the periphery (each border point links to its nearest clockwise and
#'   counterclockwise border neighbour).  Default \code{FALSE}.
#'
#' @param edge_color
#'   Character; colour of border edges.  Default \code{"gray40"}.
#'
#' @param edge_alpha
#'   Numeric; transparency of border edges (0–1).  Default \code{0.6}.
#'
#' @param edge_linewidth
#'   Numeric; line width of border edges.  Default \code{0.3}.
#'
#' @param edge_nn_k
#'   Integer; number of nearest border points to inspect when determining
#'   the clockwise/counterclockwise neighbours.  Default \code{20}.
#'
#' @param max_nn_dist
#'   Numeric; optional absolute cutoff on \code{nn_mean_dist} to remove
#'   sparse/satellite border points.  Border points with mean nearest-neighbour
#'   distance above this value are excluded.  Default \code{NULL} (no filter).
#'
#' @param nn_dist_quantile
#'   Numeric in (0, 1]; keep only border points whose \code{nn_mean_dist} is
#'   at or below this quantile of the selected border set.  Lower values
#'   remove more sparse outliers (e.g. \code{0.98} drops the sparsest 2\%).
#'   Ignored when \code{max_nn_dist} is provided.  Default \code{1} (no filter).
#'
#' @details
#'   \strong{Drawing order:}
#'   \enumerate{
#'     \item Non-border cells are drawn first as a rasterized layer
#'           (\code{ggrastr::geom_point_rast()}), coloured by cluster.
#'     \item Border cells are drawn on top as crisp vector points
#'           (\code{ggplot2::geom_point()}) in black at 3× size and
#'           full opacity.
#'     \item Cluster centroid labels are drawn last (shadowtext or repel).
#'   }
#'
#'   \strong{Colour assignment:}
#'   Uses the same shared DSatur algorithm as \code{fcs_plot_reduction()}
#'   so that neighbouring clusters receive perceptually distinct colours.
#'
#' @return
#'   If \code{return_plot = TRUE}, a ggplot2 object.
#'   If \code{return_plot = FALSE}, the plot is written to file and the
#'   function invisibly returns \code{NULL}.
#'
#' @examples
#' \dontrun{
#'   joined    <- FCSimple::fcs_join(files)
#'   reduced   <- FCSimple::fcs_reduce_dimensions(joined, algorithm = "umap")
#'   clustered <- FCSimple::fcs_cluster(reduced, algorithm = "leiden")
#'   clustered <- FCSimple::fcs_reduction_borders(clustered, "leiden")
#'
#'   # Border points from the 2D UMAP's own k-NN
#'   FCSimple::fcs_plot_borders(clustered, "leiden", "UMAP", "umap_2d")
#'
#'   # High-dimensional borders projected onto 2D UMAP
#'   FCSimple::fcs_plot_borders(
#'     clustered, "leiden", "UMAP", "high_dim",
#'     border_threshold = 0.2
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_reduce_dimensions, FCSimple::fcs_cluster,
#'   FCSimple::fcs_reduction_borders, FCSimple::fcs_plot_reduction
#'
#' @importFrom ggplot2 ggplot aes scale_color_manual labs theme_bw theme
#' @importFrom ggrastr geom_point_rast
#' @importFrom shadowtext geom_shadowtext
#' @export
fcs_plot_borders <- function(
    fcs_join_obj,
    algorithm,
    reduction,
    borders_slot,
    border_threshold = 0.3,
    point_alpha = 0.1,
    point_size = 1,
    border_size_factor = 3,
    outdir = getwd(),
    figure_width = 10,
    figure_height = 10,
    plotting_device = "pdf",
    annotate_text_size = 8,
    annotation_method = "shadowtext",
    title_size = 14,
    return_plot = FALSE,
    randomize_colors = FALSE,
    color_random_seed = 123L,
    cluster_substitute_names = NA,
    add_timestamp = TRUE,
    draw_edges = FALSE,
    edge_color = "gray40",
    edge_alpha = 0.6,
    edge_linewidth = 0.3,
    edge_nn_k = 20L,
    max_nn_dist = NULL,
    nn_dist_quantile = 1,
    border_type = c("all", "inter_cluster", "periphery"))
{
  if (!require(ggplot2, quietly = TRUE)) stop("Package 'ggplot2' is required but could not be loaded.")
  if (!require(ggrastr, quietly = TRUE)) stop("Package 'ggrastr' is required but could not be loaded.")
  if (!require(shadowtext, quietly = TRUE)) stop("Package 'shadowtext' is required but could not be loaded.")

  # -- Case-insensitive lookup ------------------------------------------------
  algorithm <- tolower(algorithm)
  borders_slot <- tolower(borders_slot)

  # -- Extract coordinates -----------------------------------------------------
  reduction_coords <- fcs_join_obj[[.resolve_reduction_slot(fcs_join_obj, reduction, 2L)]][["coordinates"]]
  if (is.null(reduction_coords))
    stop("No coordinates found under fcs_join_obj$", reduction, "$coordinates.")
  if (!is.matrix(reduction_coords))
    reduction_coords <- as.matrix(reduction_coords)
  if (ncol(reduction_coords) < 2L)
    stop("The reduction has only ", ncol(reduction_coords),
         " dimension(s); 2D is required for plotting.")

  cluster_numbers <- as.character(fcs_join_obj[[algorithm]][["clusters"]])
  if (is.null(cluster_numbers))
    stop("No cluster labels found under fcs_join_obj$", algorithm, "$clusters.")

  # -- Extract border scores --------------------------------------------------
  border_df <- fcs_join_obj[["reduction_borders"]]
  if (is.null(border_df))
    stop("No reduction_borders found in fcs_join_obj.  ",
         "Run fcs_reduction_borders() first.")
  if (!borders_slot %in% names(border_df))
    stop("Border slot '", borders_slot, "' not found in fcs_join_obj$reduction_borders.  ",
         "Available slots: ", paste(names(border_df), collapse = ", "))

  bslot <- border_df[[borders_slot]]
  if (is.null(bslot))
    stop("Slot '", borders_slot, "' is empty in reduction_borders.")

  border_data <- bslot$scores
  if (is.null(border_data) && "border_score" %in% names(bslot))
    border_data <- bslot  # legacy flat data.frame

  border_type <- match.arg(border_type)
  score_col <- switch(border_type,
    all           = "border_score",
    inter_cluster = "frac_diff",
    periphery    = "periphery_score")
  border_scores <- border_data[[score_col]]
  if (is.null(border_scores))
    stop("No '", score_col, "' column in reduction_borders$", borders_slot, ".")

  n_total <- length(cluster_numbers)

  # -- Substitute cluster names -----------------------------------------------
  if (!is.na(cluster_substitute_names[1])) {
    if (mean(names(cluster_substitute_names) %in% unique(cluster_numbers)) != 1) {
      stop("error in argument 'cluster_substitute_names': length of ",
           "'cluster_substitute_names' does not match number of clusters.")
    } else {
      message("Substituted cluster names will replace cluster numbers.")
      cluster_numbers <- as.character(cluster_substitute_names[cluster_numbers])
    }
  }

  uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]

  # -- Compute cluster centroids ----------------------------------------------
  xval <- rep(NA, length(uclus)); names(xval) <- uclus
  yval <- xval
  xmean <- xval; ymean <- xval
  for (i in seq_along(xval)) {
    cl <- names(xval)[i]
    xval[i]  <- stats::median(reduction_coords[, 1L][cluster_numbers == cl])
    yval[i]  <- stats::median(reduction_coords[, 2L][cluster_numbers == cl])
    xmean[i] <- mean(reduction_coords[, 1L][cluster_numbers == cl])
    ymean[i] <- mean(reduction_coords[, 2L][cluster_numbers == cl])
  }

  # -- Assign cluster colours (shared across plots) ---------------------------
  # Colour mapping is derived from the first-created reduction so that
  # border plots match the cluster colours of fcs_plot_reduction*().
  primary_slot <- .get_primary_reduction_slot(fcs_join_obj)
  if (randomize_colors) {
    base_colors <- fcs_color_palette(length(uclus))
    set.seed(color_random_seed)
    colorby        <- base_colors
    names(colorby) <- sample(uclus, length(uclus), replace = FALSE)
  } else {
    colorby <- .compute_cluster_color_mapping(
      fcs_join_obj, cluster_numbers, primary_slot
    )
  }

  # -- Split into border vs non-border ----------------------------------------
  is_border <- border_scores >= border_threshold

  # Optional local-density filter: remove sparse/satellite "border" points.
  # Satellite points have a high periphery_score but sit in sparse space
  # (large mean nearest-neighbour distance).  Filter by absolute cutoff
  # and/or quantile of the already-selected border points.
  if ((!is.null(max_nn_dist) || nn_dist_quantile < 1)) {
    nn_md <- border_data[["nn_mean_dist"]]
    if (is.null(nn_md) || all(is.na(nn_md))) {
      stop("This reduction_borders slot has no 'nn_mean_dist'. ",
           "Re-run fcs_reduction_borders() to compute it.")
    }
    if (!is.null(max_nn_dist)) {
      is_border <- is_border & (nn_md <= max_nn_dist)
    } else {
      cutoff <- stats::quantile(nn_md[is_border], probs = nn_dist_quantile,
                                na.rm = TRUE)
      is_border <- is_border & (nn_md <= cutoff)
    }
  }

  n_border <- sum(is_border)
  message("Border threshold ", border_threshold, ": ",
          n_border, " / ", n_total, " cells (",
          round(100 * n_border / n_total, 1), "%)")

  # Non-border data frame
  df_non <- data.frame(
    x       = reduction_coords[!is_border, 1L, drop = TRUE],
    y       = reduction_coords[!is_border, 2L, drop = TRUE],
    cluster = cluster_numbers[!is_border],
    stringsAsFactors = FALSE
  )

  # Border data frame
  df_border <- data.frame(
    x       = reduction_coords[is_border, 1L, drop = TRUE],
    y       = reduction_coords[is_border, 2L, drop = TRUE],
    cluster = cluster_numbers[is_border],
    stringsAsFactors = FALSE
  )

  # -- Build plot -------------------------------------------------------------
  coord_names <- colnames(reduction_coords)[1:2]
  if (is.null(coord_names) || any(is.na(coord_names)))
    coord_names <- c("Dim1", "Dim2")

  # Layer 1: non-border points (rasterized, coloured by cluster)
  p <- ggplot(data = df_non, mapping = aes(x = x, y = y, color = cluster)) +
    ggrastr::geom_point_rast(alpha = point_alpha, size = point_size, stroke = 0) +
    scale_color_manual(values = colorby) +
    labs(x = coord_names[1L], y = coord_names[2L])

  # Layer 2: border edges (if requested).
  # Each border point is connected to its two nearest border neighbours that
  # flank it along the periphery: the nearest clockwise and nearest
  # counterclockwise border point.  Neighbour angles around each point are
  # sorted; the largest empty angular sector marks the "outside" of the
  # cluster, and the two border points immediately flanking that gap are the
  # CW and CCW neighbours.
  if (draw_edges && n_border > 1L) {
    border_idx  <- which(is_border)
    n_border    <- length(border_idx)
    border_coords <- reduction_coords[border_idx, 1:2, drop = FALSE]

    # Nearest border points (self + up to edge_nn_k candidates)
    candidate_k  <- min(edge_nn_k + 1L, n_border)
    nn <- RANN::nn2(
      data        = border_coords,
      query       = border_coords,
      k           = candidate_k,
      treetype    = "kd",
      searchtype  = "standard"
    )
    nn_ids <- nn$nn.idx  # rows = border-local index, cols = self + candidates

    twopi <- 2 * pi
    edge_from <- integer(0L)
    edge_to   <- integer(0L)

    for (i in seq_len(n_border)) {
      cand_local <- nn_ids[i, -1L]                     # drop self
      cand_local <- cand_local[cand_local != i]        # guard against self
      if (length(cand_local) == 0L) next

      # Angles of candidate border points around the current point
      dx  <- border_coords[cand_local, 1L] - border_coords[i, 1L]
      dy  <- border_coords[cand_local, 2L] - border_coords[i, 2L]
      ang <- atan2(dy, dx)
      ord <- order(ang)
      ang <- ang[ord]
      loc <- cand_local[ord]

      # Circular gaps; largest gap = the empty sector (outside of cluster)
      gaps <- diff(c(ang, ang[1L] + twopi))
      gap_idx <- which.max(gaps)

      m <- length(loc)
      if (m == 1L) {
        # Only one candidate on this side; connect to it (its counterpart
        # is chosen by the other point's own CW/CCW logic).
        flank <- loc[1L]
        edge_from <- c(edge_from, border_idx[i])
        edge_to   <- c(edge_to,   border_idx[flank])
      } else {
        # Two flanking neighbours: one just before the gap (CCW side) and
        # one just after it (CW side).
        before <- gap_idx                       # angle just below the gap
        after  <- gap_idx %% m + 1L             # angle just above the gap
        for (flank in unique(c(before, after))) {
          edge_from <- c(edge_from, border_idx[i])
          edge_to   <- c(edge_to,   border_idx[loc[flank]])
        }
      }
    }

    # Deduplicate undirected pairs
    if (length(edge_from) > 0L) {
      key <- paste(pmin(edge_from, edge_to), pmax(edge_from, edge_to), sep = "_")
      keep <- !duplicated(key)
      edge_from <- edge_from[keep]
      edge_to   <- edge_to[keep]

      df_edges <- data.frame(
        x    = reduction_coords[edge_from, 1L],
        y    = reduction_coords[edge_from, 2L],
        xend = reduction_coords[edge_to,   1L],
        yend = reduction_coords[edge_to,   2L],
        stringsAsFactors = FALSE
      )
      p <- p + ggplot2::geom_segment(
        data      = df_edges,
        mapping   = aes(x = x, y = y, xend = xend, yend = yend),
        color     = edge_color,
        alpha     = edge_alpha,
        linewidth = edge_linewidth,
        inherit.aes = FALSE
      )
      message("Drew ", nrow(df_edges), " border edges for ", n_border,
              " border points.")
    }
  }

  # Layer 3: border points (vector, black, larger, on top)
  if (n_border > 0L) {
    p <- p + ggplot2::geom_point(
      data    = df_border,
      mapping = aes(x = x, y = y),
      color   = "black",
      alpha   = 1,
      size    = point_size * border_size_factor,
      stroke  = 0
    )
  }

  # Layer 4: centroid labels
  if (!is.na(annotate_text_size) && !is.na(annotation_method)) {
    if (annotation_method == "shadowtext") {
      centroids_df <- data.frame(
        x       = xmean,
        y       = ymean,
        cluster = names(xmean),
        stringsAsFactors = FALSE
      )
      p <- p + shadowtext::geom_shadowtext(
        data        = centroids_df,
        mapping     = aes(x = x, y = y, label = cluster),
        color       = "white",
        bg.color    = "black",
        bg.r        = 0.15,
        size        = annotate_text_size,
        fontface    = "bold",
        inherit.aes = FALSE
      )
    } else if (annotation_method == "repel") {
      if (!require(ggrepel, quietly = TRUE)) stop("Package 'ggrepel' is required but could not be loaded.")
      centroids_df <- data.frame(
        x       = xmean,
        y       = ymean,
        cluster = names(xmean),
        stringsAsFactors = FALSE
      )
      p <- p + ggrepel::geom_text_repel(
        data         = centroids_df,
        mapping      = aes(x = x, y = y, label = cluster),
        color        = "white",
        size         = annotate_text_size,
        fontface     = "bold",
        box.padding  = 0.3,
        point.padding = 0.2,
        segment.alpha = 0.4,
        inherit.aes  = FALSE
      )
    }
  }

  # Theme
  p <- p + theme_bw(base_size = title_size) +
    theme(
      legend.position = "none",
      plot.title      = element_text(face = "bold", hjust = 0.5, size = title_size)
    ) +
    ggplot2::ggtitle(label = paste0(
      reduction, " — ", algorithm, " — ", borders_slot,
      " borders (threshold = ", border_threshold, ")"))

  # -- Return or save ---------------------------------------------------------
  if (return_plot) {
    return(p)
  } else {
    if (add_timestamp) {
      fname <- paste0(
        outdir, "/", algorithm, "_", reduction, "_", borders_slot,
        "_borders_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S")
      )
    } else {
      fname <- paste0(
        outdir, "/", algorithm, "_", reduction, "_", borders_slot, "_borders"
      )
    }
    if (plotting_device == "pdf") {
      ggplot2::ggsave(
        filename = paste0(fname, ".pdf"),
        plot = p, device = "pdf", width = figure_width, height = figure_height,
        units = "in", dpi = 900, bg = "white"
      )
    } else if (plotting_device == "png") {
      ggplot2::ggsave(
        filename = paste0(fname, ".png"),
        plot = p, device = "png", width = figure_width, height = figure_height,
        units = "in", dpi = 900, bg = "white"
      )
    } else {
      stop("plotting_device must be 'pdf' or 'png'")
    }
    message("Plot saved to: ", fname, ".", plotting_device)
  }
}