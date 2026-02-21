#' @title Export Flow Cytometry Data to an FCS File
#'
#' @description
#'   Takes an FCSimple analysis object and writes its expression data (raw or
#'   transformed), with optional dimensionality‐reduction coordinates and
#'   clustering labels, into a new .fcs file. Useful for exporting combined
#'   data back into flowCore‐compatible format for downstream tools.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), optionally augmented by
#'   FCSimple::fcs_batch_correction(), FCSimple::fcs_reduce_dimensions(), and
#'   FCSimple::fcs_cluster(). Must contain at least one of:
#'   - `data`: transformed expression matrix (events × channels)
#'   - `raw`: raw (reverse‐transformed) matrix, if `data_format = "raw"` is used
#'
#' @param fcs_name
#'   Character string (without “.fcs”) to use as the base filename for the output
#'   FCS file (default `"fcs_out"`).
#'
#' @param data_format
#'   Character; which data slot to write. Options are:
#'   - `"raw"`: write `fcs_join_obj$raw` if present, otherwise falls back to
#'     `fcs_join_obj$data` (default)
#'   - `"transformed"`: write `fcs_join_obj$data`
#'
#' @param include_reductions
#'   Character vector naming reduction embeddings to include as extra columns.
#'   Each element (e.g. `"UMAP"` or `"tSNE"`) must match a named list element
#'   in `fcs_join_obj` containing `$coordinates`. Defaults to `c("UMAP","tSNE")`.
#'
#' @param include_clusterings
#'   Character vector of clustering results to include. Each element (e.g.
#'   `"leiden"`, `"flowsom"`, `"louvain"`, `"phenograph"`) must match
#'   a list element in `fcs_join_obj` with `clusters`. Cluster columns are named `<algorithm>_cluster`.
#'
#' @param subset_rows
#'   Either `"all"` (default) to write every event, or an integer vector of
#'   row indices to subset before writing.
#'
#' @param outdir
#'   Directory path where the .fcs file will be saved (default `getwd()`).
#'
#' @param include_timestamp
#'   Logical; if `TRUE` (default), appends `_YYYY-MM-DD_HHMMSS` to the filename.
#'
#' @details
#'   The function proceeds as follows:
#'   1. Validates that `outdir` exists.
#'   2. Selects the expression matrix according to `data_format`.
#'   3. Validates `subset_rows`: must be `"all"` or an integer vector with all
#'      indices in `[1, nrow(data)]`; an informative error is thrown otherwise.
#'   4. Gathers any requested embeddings from
#'      `fcs_join_obj[[ tolower(reduction) ]][["coordinates"]]`, skipping any
#'      that are missing or whose row count does not match the data matrix.
#'   5. Gathers any requested cluster assignments from
#'      `fcs_join_obj[[ tolower(algorithm) ]][["clusters"]]`, skipping any
#'      that are missing or whose length does not match the data matrix.
#'   6. Column‐binds the data, embeddings, and cluster columns into a single
#'      numeric matrix via `do.call(cbind, ...)`.
#'   7. Constructs a `flowFrame` and, if `subset_rows != "all"`, subsets it.
#'   8. Writes to disk via `flowCore::write.FCS()` using `file.path()` for
#'      cross‐platform path construction.
#'
#' @return
#'   Invisibly returns `NULL`. The main effect is the creation of an .fcs file
#'   at `file.path(outdir, paste0(fcs_name, [optional _timestamp], ".fcs"))`.
#'
#' @examples
#' \dontrun{
#'   # Export transformed data only
#'   joined <- FCSimple::fcs_join(file_paths)
#'   FCSimple::fcs_write.FCS(
#'     joined,
#'     fcs_name = "my_analysis",
#'     data_format = "transformed",
#'     include_reductions = c("UMAP"),
#'     include_clusterings = c("leiden"),
#'     outdir = "~/results"
#'   )
#'
#'   # Export raw values + all reductions and clusters, subset first 1000 cells
#'   FCSimple::fcs_write.FCS(
#'     joined,
#'     data_format = "raw",
#'     subset_rows = 1:1000,
#'     include_timestamp = FALSE
#'   )
#' }
#'
#' @seealso
#'   flowCore::write.FCS, FCSimple::fcs_join,
#'   FCSimple::fcs_reduce_dimensions, FCSimple::fcs_cluster
#'
#' @importFrom flowCore write.FCS
#' @export
fcs_write.FCS <- function(fcs_join_obj,
                          fcs_name = "fcs_out",
                          data_format = c("raw","transformed"),
                          include_reductions = c("UMAP","tSNE"),
                          include_clusterings = c("leiden","flowsom","louvain","phenograph"),
                          subset_rows = "all", # either "all" or a numeric vector of positions to write to file
                          outdir = getwd(),
                          include_timestamp = TRUE)
{
  require(flowCore)

  # --- outdir validation ---
  if (!dir.exists(outdir)) {
    stop(paste0("'outdir' does not exist: ", outdir))
  }

  # --- select expression matrix ---
  if (length(data_format) == 2) {
    data_format <- "raw"
  }
  if (data_format == "raw") {
    if (!"raw" %in% names(fcs_join_obj)) {
      warning("'raw' format not found in object. Consider adding a 'raw' element with reverse-transformed values. Using transformed values stored in $data.")
      data_incl <- fcs_join_obj[["data"]]
    } else {
      data_incl <- fcs_join_obj[["raw"]]
    }
  } else if (data_format == "transformed") {
    data_incl <- fcs_join_obj[["data"]]
  } else {
    stop("error in argument 'data_format': use either 'raw' or 'transformed'")
  }
  if (is.null(data_incl)) {
    stop("Expression matrix is NULL. Verify that the input object has a 'data' (or 'raw') entry.")
  }
  n_rows <- nrow(data_incl)

  # --- validate subset_rows ---
  if (!identical(subset_rows, "all")) {
    if (!is.numeric(subset_rows) && !is.integer(subset_rows)) {
      stop("'subset_rows' must be \"all\" or an integer/numeric vector of row indices.")
    }
    out_of_range <- subset_rows[subset_rows < 1 | subset_rows > n_rows]
    if (length(out_of_range) > 0) {
      stop(paste0("'subset_rows' contains ", length(out_of_range),
                  " out-of-range index/indices (max valid: ", n_rows, ")."))
    }
  }

  fcs_include <- list(data = data_incl)

  # --- reductions ---
  if (length(include_reductions) != 0) {
    reduction_parts <- list()
    for (r in include_reductions) {
      key <- tolower(r)
      if (key %in% names(fcs_join_obj)) {
        coords <- as.matrix(fcs_join_obj[[key]][["coordinates"]])
        if (nrow(coords) != n_rows) {
          warning(paste0("Reduction '", r, "' has ", nrow(coords),
                         " rows but data has ", n_rows, " rows; skipping."))
        } else {
          reduction_parts[[length(reduction_parts) + 1]] <- coords
        }
      } else {
        warning(paste0("Reduction '", r, "' not found in object and will not be included."))
      }
    }
    if (length(reduction_parts) > 0) {
      fcs_include[["reductions"]] <- do.call(cbind, reduction_parts)
    }
  }

  # --- clusterings ---
  if (length(include_clusterings) != 0) {
    cluster_parts <- list()
    cluster_names <- character(0)
    for (a in include_clusterings) {
      key <- tolower(a)
      if (key %in% names(fcs_join_obj)) {
        clus <- fcs_join_obj[[key]][["clusters"]]
        if (length(clus) != n_rows) {
          warning(paste0("Clustering '", a, "' has ", length(clus),
                         " elements but data has ", n_rows, " rows; skipping."))
        } else {
          cluster_parts[[length(cluster_parts) + 1]] <- as.numeric(clus)
          cluster_names <- c(cluster_names, paste0(key, "_cluster"))
        }
      } else {
        warning(paste0("Clustering '", a, "' not found in object and will not be included."))
      }
    }
    if (length(cluster_parts) > 0) {
      cluster_mat <- matrix(unlist(cluster_parts), ncol = length(cluster_parts))
      colnames(cluster_mat) <- cluster_names
      fcs_include[["clusterings"]] <- cluster_mat
    }
  }

  # --- combine ---
  final_return <- as.matrix(do.call(cbind, fcs_include))

  # --- build flowFrame ---
  out_ff <- new("flowFrame", exprs = final_return)
  if (!identical(subset_rows, "all")) {
    out_ff <- out_ff[subset_rows, ]
  }

  # --- write ---
  base_name <- gsub("\\.fcs$", "", fcs_name)
  if (include_timestamp) {
    out_file <- file.path(outdir, paste0(base_name, "_", strftime(Sys.time(), "%Y-%m-%d_%H%M%S"), ".fcs"))
  } else {
    out_file <- file.path(outdir, paste0(base_name, ".fcs"))
  }
  flowCore::write.FCS(x = out_ff, filename = out_file)
  invisible(NULL)
}
