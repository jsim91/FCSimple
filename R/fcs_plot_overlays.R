#' @title Project Marker Expression onto a 2D Embedding
#'
#' @description
#'   Creates a multi‐panel PDF of channel (parameter) expression projected onto
#'   a 2D reduction embedding (UMAP or tSNE). For each specified parameter,
#'   events are optionally trimmed of outliers, subsampled, and plotted as a
#'   rasterized scatter colored by expression intensity. Panels are arranged
#'   in 2×2 grids and written to disk.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), optionally augmented by
#'   fcs_batch_correction() and fcs_reduce_dimensions(), containing:
#'   - `data`: raw or transformed expression matrix (events × channels)
#'   - `batch_correction$data` (if present)
#'   - a `umap` or `tsne` element with `$coordinates` (events × 2)
#'
#' @param override_correction
#'   Logical; if `TRUE` (default), use `fcs_join_obj$data` even when
#'   batch correction exists. If `FALSE`, uses
#'   `fcs_join_obj$batch_correction$data`.
#'
#' @param reduction
#'   Character; embedding to use. Must be either `"UMAP"` (default) or `"tSNE"`.
#'
#' @param parameters
#'   Character vector of channel names to project. If `"all"` (default),
#'   all columns in the selected data matrix are plotted.
#'
#' @param outdir
#'   Character; path to an existing directory where the PDF will be saved.
#'   Defaults to `getwd()`.
#'
#' @param sample_size
#'   Integer; maximum number of events to sample before plotting. Sampling is
#'   performed **once** on the common set of cells so that all parameter panels
#'   display exactly the same events, enabling direct cross‐marker comparison
#'   (default 50000).
#'
#' @param point_size
#'   Numeric; point size for the rasterized scatter (`geom_point_rast`)
#'   (default 0.8).
#'
#' @param trim_outliers
#'   Logical; if `TRUE` (default), expression values beyond
#'   `trim_quantile` are winsorized (clamped) to the quantile boundaries
#'   rather than removing rows, so that every panel retains the identical
#'   set of cells.
#'
#' @param trim_quantile
#'   Numeric; quantile threshold for outlier winsorization (default 0.01).
#'
#' @param point_alpha
#'   Numeric; point transparency (alpha) for scatter (default 0.1).
#'
#' @param force_xlim
#'   Numeric vector of length 2 to fix the x‐axis limits, or `FALSE`
#'   (default) to use the data range.
#'
#' @param force_ylim
#'   Numeric vector of length 2 to fix the y‐axis limits, or `FALSE`
#'   (default) to use the data range.
#'
#' @details
#'   1. Determines which expression matrix to use based on
#'      `override_correction`.
#'   2. Extracts UMAP or tSNE coordinates from
#'      `fcs_join_obj[[ tolower(reduction) ]][["coordinates"]]`.
#'   3. Selects the requested `parameters` (all by default).
#'   4. **Downsamples once** to a common set of up to `sample_size` cells so
#'      that every parameter panel shows the identical cell population.
#'   5. For each parameter:
#'      - Optionally winsorizes expression outliers at the specified quantiles
#'        (clamps values; no cells are removed).
#'      - Creates a rasterized scatter plot colored by expression
#'        (`scale_color_viridis(option="D")`).
#'      - Honors any fixed `force_xlim`/`force_ylim`.
#'   6. Arranges plots into pages of 2×2 panels via ggpubr::ggarrange().
#'   7. Saves a timestamped PDF named
#'      `<reduction>_parameter_projections_<YYYY-MM-DD_HHMMSS>.pdf`.
#'
#' @return
#'   Invisibly returns `NULL`. The side effect is a multi‐page PDF written
#'   to `outdir`.
#'
#' @examples
#' \dontrun{
#'   # Assume joined, reduced and batch-corrected object
#'   joined <- FCSimple::fcs_join(files)
#'   reduced <- FCSimple::fcs_reduce_dimensions(joined, algorithm = "umap")
#'   corrected <- FCSimple::fcs_batch_correction(reduced)
#'
#'   # Project all channels onto UMAP
#'   FCSimple::fcs_plot_overlays(
#'     corrected,
#'     reduction = "UMAP",
#'     outdir    = "~/results"
#'   )
#'
#'   # Project only CD3 and CD19
#'   FCSimple::fcs_plot_overlays(
#'     corrected,
#'     parameters = c("CD3","CD19"),
#'     sample_size = 30000
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_reduce_dimensions, FCSimple::fcs_plot_reduction,
#'   ggplot2::ggplot, viridis::scale_color_viridis,
#'   ggrastr::geom_point_rast, ggpubr::ggarrange
#'
#' @importFrom ggplot2 ggplot aes theme_bw theme xlim ylim ggtitle element_text
#'   element_blank unit margin guide_colorbar geom_blank
#' @importFrom viridis scale_color_viridis
#' @importFrom ggrastr geom_point_rast
#' @importFrom ggpubr ggarrange
#' @importFrom gridExtra marrangeGrob
#' @importFrom scales squish
#' @export
fcs_plot_overlays <- function(fcs_join_obj,
                                   override_correction = TRUE,
                                   reduction = c("UMAP","tSNE"),
                                   parameters = "all",
                                   outdir = getwd(),
                                   sample_size = 50000,
                                   point_size = 0.8,
                                   point_alpha = 0.1,
                                   trim_outliers = TRUE,
                                   trim_quantile = 0.01,
                                   force_xlim = FALSE,
                                   force_ylim = FALSE)
{
  # ---- input validation & data extraction -----------------------------------
  if (length(reduction) != 1L) {
    stop("error in argument 'reduction': use either 'UMAP' or 'tSNE'")
  }

  if ("batch_correction" %in% names(fcs_join_obj)) {
    if (override_correction) {
      message("batch_correction found in fcs_join_obj and override_correction set to TRUE. ",
              "Using fcs_join_obj[[\"data\"]] for projections. ",
              "To use batch-corrected features, set 'override_correction' to FALSE.")
      join_data <- fcs_join_obj[["data"]]
    } else {
      message("batch_correction found in fcs_join_obj and override_correction set to FALSE. ",
              "Using fcs_join_obj[[\"batch_correction\"]][[\"data\"]] for projections. ",
              "To use original features, set 'override_correction' to TRUE.")
      join_data <- fcs_join_obj[["batch_correction"]][["data"]]
    }
  } else {
    join_data <- fcs_join_obj[["data"]]
  }

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  if (is.null(reduction_coords) || nrow(reduction_coords) < 1L) {
    stop("error in 'reduction': unable to find specified reduction")
  }

  if (nrow(join_data) != nrow(reduction_coords)) {
    stop("Mismatch: expression data has ", nrow(join_data), " rows but reduction has ",
         nrow(reduction_coords), " rows.")
  }

  # ---- parameter selection --------------------------------------------------
  if (identical(parameters, "all")) {
    include_params <- colnames(join_data)
  } else {
    target_params <- intersect(parameters, colnames(join_data))
    if (length(target_params) == 0L) {
      stop("error in argument 'parameters': none of the parameters requested were found")
    }
    include_params <- target_params
  }

  # ---- determine reduction column names once --------------------------------
  red_lower <- tolower(reduction)
  if (grepl("umap", red_lower, fixed = TRUE)) {
    red_x <- "UMAP1"; red_y <- "UMAP2"
  } else {
    red_x <- "tSNE1"; red_y <- "tSNE2"
  }

  # ---- build a single data.frame & downsample ONCE --------------------------
  #  Combine reduction coordinates with all requested expression columns so
  #  that every parameter panel draws from exactly the same set of cells.
  plot_df <- as.data.frame(cbind(reduction_coords, join_data[, include_params, drop = FALSE]))
  colnames(plot_df)[1:2] <- c(red_x, red_y)

  total_n <- nrow(plot_df)
  if (total_n > sample_size) {
    set.seed(123L)
    keep_idx <- sample(seq_len(total_n), size = sample_size, replace = FALSE)
    plot_df <- plot_df[keep_idx, , drop = FALSE]
  }

  # ---- internal: build a single-panel plot ----------------------------------
  build_panel <- function(df, param_name,
                          pts = point_size, palpha = point_alpha,
                          tr_out = trim_outliers, tr_q = trim_quantile) {
    # Winsorize expression values to the quantile range instead of removing
    # rows.  This keeps the cell set identical across all panels.
    xvals <- df[[param_name]]
    if (isTRUE(tr_out)) {
      qlim <- quantile(xvals, probs = c(tr_q, 1 - tr_q), na.rm = TRUE)
      xvals <- pmax(pmin(xvals, qlim[2]), qlim[1])
    }

    plt <- ggplot(data = df, mapping = aes(x = .data[[red_x]],
                                           y = .data[[red_y]],
                                           color = xvals)) +
      geom_point_rast(size = pts, pch = 19, alpha = palpha) +
      scale_color_viridis(
        option   = "D",
        oob      = scales::squish,               # safety-net for any remaining OOB values
        guide    = guide_colorbar(
          title.position = "top",
          frame.colour   = "black",
          ticks.colour   = "black",
          frame.linewidth = 0.4,
          draw.ulim      = TRUE,
          draw.llim      = TRUE,
          label.theme    = element_text(angle = 90, vjust = 0.5, size = 16),
          position       = "bottom",
          barwidth       = 10,
          barheight      = 1
        )
      ) +
      ggtitle(param_name) +
      theme_bw(base_size = 18) +
      theme(
        axis.title     = element_blank(),
        axis.text      = element_blank(),
        axis.ticks     = element_blank(),
        legend.title   = element_blank(),
        plot.title     = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.spacing  = unit(0, "pt"),
        legend.margin   = margin(t = -5, r = 0, b = 20, l = 0),
        plot.margin     = margin(t = 5,  r = 5, b = 0,  l = 5)
      )

    if (is.numeric(force_xlim)) plt <- plt + xlim(force_xlim)
    if (is.numeric(force_ylim)) plt <- plt + ylim(force_ylim)

    plt
  }

  # ---- generate one plot per parameter (same cells in every panel) ----------
  intens_plots <- lapply(include_params, function(p) build_panel(plot_df, p))

  # ---- arrange into 2×2 pages -----------------------------------------------
  #  Create a true blank placeholder for pages that don't have 4 panels.
  void_plot <- ggplot() + geom_blank() + theme_void()

  n_plots   <- length(intens_plots)
  n_pages   <- ceiling(n_plots / 4)
  n_padded  <- n_pages * 4L

  # Pad with void plots so every page has exactly 4 panels.
  if (n_padded > n_plots) {
    intens_plots[seq.int(n_plots + 1L, n_padded)] <- list(void_plot)
  }

  arranged_list <- vector("list", length = n_pages)
  for (i in seq_len(n_pages)) {
    idx <- seq.int((i - 1L) * 4L + 1L, i * 4L)
    arranged_list[[i]] <- ggpubr::ggarrange(
      plotlist = intens_plots[idx],
      nrow = 2, ncol = 2
    )
  }

  # ---- save PDF -------------------------------------------------------------
  reduc_label <- if (grepl("umap", red_lower, fixed = TRUE)) "UMAP" else "tSNE"
  out_file <- paste0(reduc_label, "_parameter_projections_",
                     strftime(Sys.time(), "%Y-%m-%d_%H%M%S"), ".pdf")

  ggsave(
    filename = out_file,
    plot     = gridExtra::marrangeGrob(grobs = arranged_list, nrow = 1, ncol = 1, top = ""),
    device   = "pdf",
    path     = outdir,
    width    = 12,
    height   = 12,
    units    = "in",
    dpi      = 900
  )
}