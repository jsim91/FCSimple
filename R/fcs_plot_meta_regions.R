#' @title Visualize Meta-Region Integrity on a 2D Embedding
#'
#' @description
#'   Colours a 2D reduction embedding by each cell's meta-region, with optional
#'   overlays highlighting meta-regions whose composition is poorly preserved
#'   (low retention / high mixing) after projection from high dimensions.
#'
#' @param fcs_join_obj
#'   A list produced by \code{FCSimple::fcs_embedding_neighbors()} and
#'   \code{FCSimple::fcs_meta_regions()}.
#'
#' @param reduction
#'   Character; the 2D reduction to plot (e.g. \code{"umap"}).
#'   Case-insensitive.
#'
#' @param metric
#'   Character; per-meta-region integrity metric to highlight:
#'   \code{"retention"} (default) or \code{"mixing"}.
#'
#' @param highlight_top
#'   Integer; number of worst meta-regions (lowest retention or highest mixing)
#'   to highlight as enlarged markers.  Default \code{20}.
#'
#' @param point_alpha
#'   Numeric; point transparency for background cells.  Default \code{0.2}.
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
#' @details
#'   Background cells are coloured by meta-region id (categorical palette).
#'   Then the \code{highlight_top} worst meta-regions according to
#'   \code{metric} are drawn as enlarged black-outlined points along with their
#'   meta-region ids, so the viewer can locate where high-dimensional coherent
#'   patches break apart in the low-dimensional layout.
#'
#' @return
#'   If \code{return_plot = TRUE}, a ggplot object; otherwise the plot is
#'   written to \code{outdir} and \code{NULL} is returned invisibly.
#'
#' @examples
#' \dontrun{
#'   obj <- fcs_meta_regions(obj, resolution_parameter = 5, n_hops = 2)
#'   fcs_plot_meta_regions(obj, "umap", metric = "retention")
#' }
#'
#' @seealso
#'   FCSimple::fcs_meta_regions, FCSimple::fcs_embedding_neighbors
#'
#' @importFrom ggplot2 ggplot aes scale_color_manual labs theme_bw theme
#'   ggsave geom_point
#' @importFrom ggrastr geom_point_rast
#' @export
fcs_plot_meta_regions <- function(
    fcs_join_obj,
    reduction,
    metric = "retention",
    highlight_top = 20,
    point_alpha = 0.2,
    outdir = getwd(),
    figure_width = 10,
    figure_height = 10,
    plotting_device = "pdf",
    add_timestamp = TRUE,
    return_plot = FALSE)
{
  if (!require(ggplot2, quietly = TRUE)) stop("Package 'ggplot2' is required but could not be loaded.")
  if (!require(ggrastr, quietly = TRUE)) stop("Package 'ggrastr' is required but could not be loaded.")

  metric <- match.arg(metric, c("retention", "mixing"))
  mr <- fcs_join_obj$meta_regions
  if (is.null(mr))
    stop("No '$meta_regions' found. Run fcs_meta_regions() first.")

  red <- tolower(reduction)
  slot <- paste0(red, "_2d")
  if (!slot %in% names(mr$regions))
    slot <- grep(paste0("^", red, "_2d$"), names(mr$regions), value = TRUE)[1L]
  if (is.na(slot) || !slot %in% names(mr$regions))
    stop("No meta-region data for reduction '", red, "' (2D).")

  # Coordinates
  if (slot %in% names(fcs_join_obj)) {
    coords <- fcs_join_obj[[slot]]$coordinates
  } else if (red %in% names(fcs_join_obj)) {
    coords <- fcs_join_obj[[red]]$coordinates
  } else {
    stop("No coordinates for reduction '", red, "' found.")
  }
  if (!is.matrix(coords)) coords <- as.matrix(coords)

  membership <- mr$membership
  region_stat <- mr$regions[[slot]]

  df <- data.frame(
    x = coords[, 1L], y = coords[, 2L],
    meta = factor(membership),
    stringsAsFactors = FALSE
  )

  p <- ggplot(df, aes(x = x, y = y, color = meta)) +
    ggrastr::geom_point_rast(alpha = point_alpha, size = 0.5, stroke = 0) +
    scale_color_manual(values = fcs_color_palette(length(levels(df$meta)))) +
    labs(x = colnames(coords)[1L], y = colnames(coords)[2L]) +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))

  # Highlight worst meta-regions
  if (highlight_top > 0 && highlight_top < nrow(region_stat)) {
    worst <- if (metric == "retention")
      head(region_stat[order(region_stat$retention), ], highlight_top)
    else
      head(region_stat[order(-region_stat$mixing), ], highlight_top)

    worst_ids <- worst$meta_region
    keep <- membership %in% worst_ids
    hdf <- data.frame(
      x = coords[keep, 1L], y = coords[keep, 2L],
      meta = factor(membership[keep]),
      stringsAsFactors = FALSE
    )
    p <- p + ggplot2::geom_point(
      data = hdf, mapping = aes(x = x, y = y, color = meta),
      size = 1.5, stroke = 0.3, alpha = 1)
  }

  p <- p + ggplot2::ggtitle(paste0(
    reduction, " — meta-region ", metric, " (top ", highlight_top, " highlighted)"))

  if (return_plot) return(p)

  if (add_timestamp) {
    fname <- paste0(outdir, "/", slot, "_meta_regions_", metric, "_",
                    strftime(Sys.time(), "%Y-%m-%d_%H%M%S"))
  } else {
    fname <- paste0(outdir, "/", slot, "_meta_regions_", metric)
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