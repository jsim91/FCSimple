#' @title Interactive 3D Reduction Embedding Plot
#'
#' @description
#'   Creates an interactive 3D scatter plot of a dimensionality reduction
#'   embedding (UMAP or t-SNE) with points coloured by cluster membership.
#'   Uses plotly for WebGL-accelerated rendering with rotation, zoom, and
#'   hover interactivity.  Cluster labels are placed at 3D centroids and
#'   remain readable as the user rotates the view.  For large datasets,
#'   automatic stratified downsampling preserves rare populations while
#'   keeping the plot responsive.
#'
#' @param fcs_join_obj
#'   A list containing reduction coordinates and clustering results, as
#'   produced by FCSimple::fcs_reduce_dimensions() and
#'   FCSimple::fcs_cluster().  The reduction must have been run with
#'   `n_components = 3` (or higher).
#'
#' @param algorithm
#'   Character; clustering result to colour by (e.g. `"leiden"`, `"flowsom"`).
#'
#' @param reduction
#'   Character; dimensionality-reduction to plot. Either `"UMAP"` or `"tSNE"`.
#'
#' @param point_alpha
#'   Numeric; point opacity (0–1).  Default `0.5`.
#'
#' @param point_size
#'   Numeric; point marker size.  Default `3`.
#'
#' @param max_points
#'   Integer; maximum number of points to render.  When the dataset exceeds
#'   this limit, cells are downsampled proportionally within each cluster
#'   (stratified sampling) so that rare populations are preserved.
#'   Default `50000`.
#'
#' @param downsample_seed
#'   Integer; seed for reproducible downsampling.  Default `123`.
#'
#' @param annotate_text_size
#'   Numeric; font size for centroid cluster labels.  Set `NA` to hide
#'   labels.  Default `10`.
#'
#' @param randomize_colors
#'   Logical; if `TRUE`, shuffles cluster-colour assignment (default `FALSE`).
#'
#' @param color_random_seed
#'   Integer; seed for random colour assignment when `randomize_colors = TRUE`
#'   (default `123`).
#'
#' @param color_clusters
#'   Logical; if `TRUE` (default), colour points by cluster; if `FALSE`, all
#'   points are drawn in a single colour.
#'
#' @param cluster_substitute_names
#'   Named character vector; optional mapping from original cluster IDs to
#'   replacement labels.  The names of this vector must exactly match the
#'   unique cluster IDs in
#'   `fcs_join_obj[[tolower(algorithm)]][["clusters"]]`; a mismatch will
#'   trigger an error.  Default is `NA` (no substitution).
#'
#' @param return_plot
#'   Logical; if `TRUE` (default), returns the plotly widget object (which
#'   renders interactively in RStudio, R Markdown, or Quarto).  If `FALSE`,
#'   saves a self-contained `.html` file to `outdir` and invisibly returns
#'   the file path.
#'
#' @param outdir
#'   Character; directory path to save the `.html` file when
#'   `return_plot = FALSE`.  Default `getwd()`.
#'
#' @param add_timestamp
#'   Logical; if `TRUE` (default), appends a timestamp
#'   (`_YYYY-MM-DD_HHMMSS`) to the output filename.
#'
#' @param figure_width
#'   Numeric; width in pixels for the plotly widget.  Default `900`.
#'
#' @param figure_height
#'   Numeric; height in pixels for the plotly widget.  Default `700`.
#'
#' @details
#'   \enumerate{
#'     \item Extracts the 3D embedding coordinates and cluster IDs.
#'     \item Computes 3D cluster centroids (median for colour optimisation,
#'           mean for label placement).
#'     \item Assigns colours using the same space-aware DSatur algorithm as
#'           `fcs_plot_reduction()`: a k-NN graph is built on 3D centroid
#'           distances, and colours are greedily assigned so that
#'           neighbouring clusters receive perceptually distinct colours.
#'     \item If the number of cells exceeds `max_points`, performs
#'           stratified downsampling proportional to cluster size.
#'     \item Renders an interactive plotly 3D scatter plot with per-cluster
#'           colouring, centroid text labels, and hover tooltips.
#'   }
#'
#'   **Output behaviour:**
#'   \itemize{
#'     \item `return_plot = TRUE` (default): returns the plotly widget.
#'           Displays interactively in RStudio Viewer, knitr/Quarto HTML
#'           output, or the R console.
#'     \item `return_plot = FALSE`: writes a self-contained `.html` file
#'           to `outdir` using `htmlwidgets::saveWidget()`.  The file can
#'           be opened in any modern browser with no R dependency.
#'   }
#'
#' @return
#'   If `return_plot = TRUE`, a plotly widget object.
#'   If `return_plot = FALSE`, invisibly returns the full path to the
#'   saved `.html` file.
#'
#' @examples
#' \dontrun{
#'   joined   <- FCSimple::fcs_join(files)
#'   reduced  <- FCSimple::fcs_reduce_dimensions(
#'     joined, algorithm = "umap", n_components = 3
#'   )
#'   clustered <- FCSimple::fcs_cluster(reduced, algorithm = "leiden")
#'
#'   # Interactive 3D plot in RStudio Viewer
#'   fcs_plot_reduction_3d(clustered, "leiden", "UMAP")
#'
#'   # Save as standalone HTML
#'   fcs_plot_reduction_3d(
#'     clustered, "leiden", "UMAP",
#'     return_plot = FALSE,
#'     max_points = 30000
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_reduce_dimensions, FCSimple::fcs_cluster,
#'   FCSimple::fcs_plot_reduction, plotly::plot_ly
#'
#' @importFrom plotly plot_ly layout add_trace
#' @importFrom htmlwidgets saveWidget
#' @export
fcs_plot_reduction_3d <- function(
    fcs_join_obj,
    algorithm,
    reduction,
    point_alpha = 0.5,
    point_size = 3,
    max_points = 50000L,
    downsample_seed = 123L,
    annotate_text_size = 14,
    randomize_colors = FALSE,
    color_random_seed = 123L,
    color_clusters = TRUE,
    cluster_substitute_names = NA,
    return_plot = TRUE,
    outdir = getwd(),
    add_timestamp = TRUE,
    figure_width = 1100,
    figure_height = 850)
{
  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required for fcs_plot_reduction_3d(). ",
         "Install with: install.packages('plotly')")
  if (!requireNamespace("htmlwidgets", quietly = TRUE))
    stop("Package 'htmlwidgets' is required for saving. ",
         "Install with: install.packages('htmlwidgets')")

  # -- Extract coordinates and clusters ---------------------------------------
  reduction_coords <- fcs_join_obj[[.resolve_reduction_slot(fcs_join_obj, reduction, 3L)]][["coordinates"]]
  n_dims <- ncol(reduction_coords)
  if (n_dims < 3L)
    stop("The reduction '", reduction, "' has only ", n_dims,
         " dimension(s).  Run fcs_reduce_dimensions() with n_components >= 3 ",
         "to use fcs_plot_reduction_3d().")
  if (n_dims > 3L)
    warning("The reduction has ", n_dims,
            " dimensions; only the first 3 will be plotted.")

  cluster_numbers <- as.character(
    fcs_join_obj[[tolower(algorithm)]][["clusters"]]
  )

  # -- Substitute cluster names -----------------------------------------------
  if (!is.na(cluster_substitute_names[1])) {
    if (mean(names(cluster_substitute_names) %in% unique(cluster_numbers)) != 1) {
      stop("error in argument 'cluster_substitute_names': length of ",
           "'cluster_substitute_names' does not match number of clusters.")
    } else {
      message("Substituted cluster names will replace cluster numbers. ",
              "You may use '\\n' character to create multi-line annotations.")
      cluster_numbers <- as.character(
        cluster_substitute_names[cluster_numbers]
      )
    }
  }

  uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  n_clusters <- length(uclus)

  # -- Assign cluster colours (shared across 2D and 3D) -----------------------
  # Colour mapping is derived from the first-created reduction so that plots
  # of the same clusters are consistent across embeddings.
  primary_slot <- .get_primary_reduction_slot(fcs_join_obj)
  if (randomize_colors) {
    base_colors <- fcs_color_palette(n_clusters)
    set.seed(color_random_seed)
    colorby        <- base_colors
    names(colorby) <- sample(uclus, n_clusters, replace = FALSE)
  } else {
    colorby <- .compute_cluster_color_mapping(
      fcs_join_obj, cluster_numbers, primary_slot
    )
  }

  # -- Stratified downsampling ------------------------------------------------
  n_total <- length(cluster_numbers)
  keep_idx <- seq_len(n_total)
  if (n_total > max_points) {
    message("Downsampling ", n_total, " -> ", max_points,
            " points (stratified by cluster)")
    set.seed(downsample_seed)
    cluster_counts <- table(cluster_numbers)
    cluster_frac   <- as.numeric(cluster_counts) / n_total
    alloc          <- pmax(round(max_points * cluster_frac), 1L)
    # Trim overallocation
    while (sum(alloc) > max_points) {
      over <- which.max(alloc)
      alloc[over] <- alloc[over] - 1L
    }
    names(alloc) <- names(cluster_counts)

    keep_idx <- unlist(lapply(names(alloc), function(cl) {
      idx <- which(cluster_numbers == cl)
      if (length(idx) <= alloc[cl]) return(idx)
      sample(idx, size = alloc[cl], replace = FALSE)
    }), use.names = FALSE)
    keep_idx <- sort(keep_idx)

    reduction_coords <- reduction_coords[keep_idx, , drop = FALSE]
    cluster_numbers   <- cluster_numbers[keep_idx]
  }

  # -- Compute 3D centroid means for label placement --------------------------
  centroid_mean <- matrix(NA_real_, nrow = n_clusters, ncol = 3L)
  rownames(centroid_mean) <- uclus

  for (i in seq_len(n_clusters)) {
    cl_name <- uclus[i]
    idx     <- which(cluster_numbers == cl_name)
    centroid_mean[i, ] <- colMeans(
      reduction_coords[idx, 1:3, drop = FALSE]
    )
  }

  # -- Build plotly 3D scatter ------------------------------------------------
  coord_names <- colnames(reduction_coords)[1:3]
  if (is.null(coord_names) || any(is.na(coord_names)))
    coord_names <- c("Dim1", "Dim2", "Dim3")

  df <- data.frame(
    x       = reduction_coords[, 1L],
    y       = reduction_coords[, 2L],
    z       = reduction_coords[, 3L],
    cluster = cluster_numbers,
    stringsAsFactors = FALSE
  )

  if (color_clusters) {
    p <- plotly::plot_ly(
      df,
      x = ~x, y = ~y, z = ~z,
      color = ~cluster,
      colors = unname(colorby),
      text = ~cluster,
      type = "scatter3d",
      mode = "markers",
      width = figure_width,
      height = figure_height,
      marker = list(
        size = point_size,
        opacity = point_alpha
      ),
      hovertemplate = paste0(
        "<b>%{text}</b><br>",
        coord_names[1L], ": %{x:.3f}<br>",
        coord_names[2L], ": %{y:.3f}<br>",
        coord_names[3L], ": %{z:.3f}<extra></extra>"
      )
    )
  } else {
    p <- plotly::plot_ly(
      df,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d",
      mode = "markers",
      width = figure_width,
      height = figure_height,
      marker = list(
        size = point_size,
        opacity = point_alpha,
        color = "#1f77b4"
      ),
      hovertemplate = paste0(
        coord_names[1L], ": %{x:.3f}<br>",
        coord_names[2L], ": %{y:.3f}<br>",
        coord_names[3L], ": %{z:.3f}<extra></extra>"
      )
    )
  }

  # -- Add centroid labels (single black trace, one text label per cluster) ---
  # Labels are drawn in black with a heavy font for maximum readability.  A
  # single text trace avoids the WebGL z-ordering issues that come with
  # multi-trace shadow/outline workarounds in plotly.
  if (!is.na(annotate_text_size) && n_clusters > 0L) {
    lbl_df <- data.frame(
      x     = centroid_mean[, 1L],
      y     = centroid_mean[, 2L],
      z     = centroid_mean[, 3L],
      label = rownames(centroid_mean),
      stringsAsFactors = FALSE
    )
    p <- plotly::add_trace(
      p,
      data        = lbl_df,
      x           = ~x,
      y           = ~y,
      z           = ~z,
      type        = "scatter3d",
      mode        = "text",
      text        = ~label,
      textfont    = list(
        size   = annotate_text_size,
        color  = "#000000",
        family = "Arial Black"
      ),
      hoverinfo   = "none",
      showlegend  = FALSE,
      inherit     = FALSE
    )
  }

  # -- Layout -----------------------------------------------------------------
  p <- plotly::layout(
    p,
    title = list(
      text = paste0(reduction, " (", algorithm, ")"),
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = coord_names[1L]),
      yaxis = list(title = coord_names[2L]),
      zaxis = list(title = coord_names[3L]),
      camera = list(
        eye = list(x = 1.5, y = 1.5, z = 1.2)
      )
    ),
    legend = list(
      title = list(text = "Cluster"),
      itemsizing = "constant"
    )
  )

  # -- Return or save ---------------------------------------------------------
  if (return_plot) {
    return(p)
  } else {
    if (add_timestamp) {
      fname <- paste0(
        outdir, "/", tolower(algorithm), "_", tolower(reduction),
        "_3d_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S"), ".html"
      )
    } else {
      fname <- paste0(
        outdir, "/", tolower(algorithm), "_", tolower(reduction),
        "_3d.html"
      )
    }
    htmlwidgets::saveWidget(
      widget   = p,
      file     = fname,
      title    = paste0(algorithm, " ", reduction, " 3D"),
      selfcontained = TRUE
    )
    message("Plot saved to: ", fname)
    invisible(fname)
  }
}