#' @title Prepare FCSimple Object for FCView App Upload
#'
#' @description
#'   Prepares a flow cytometry analysis object for upload to the FCView Shiny
#'   application by removing non‐essential fields, optionally downsampling,
#'   and standardizing the clustering algorithm name to 'cluster'. Can save
#'   the prepared object as an .RData file.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple workflow functions (e.g., fcs_join,
#'   fcs_cluster, fcs_reduce_dimensions), containing at minimum:
#'   - `data`: numeric matrix of events × channels
#'   - `source`: character vector of sample identifiers
#'   - `metadata`: data frame with sample‐level metadata
#'   - At least one clustering algorithm result (e.g., `leiden`, `flowsom`)
#'     whose slot contains pre-computed `frequency` and `counts` tables
#'     (produced by `FCSimple::fcs_calculate_abundance()`). Both tables
#'     **must** be present; the function will stop with an error if either
#'     is missing. They are carried through as-is and never recalculated.
#'
#' @param downsample_size
#'   Integer or `NULL`; maximum number of cells to retain. If `NULL` (default),
#'   no downsampling is performed. If specified, randomly samples cells and
#'   updates all associated fields (data, source, run_date, cluster assignments,
#'   UMAP/tSNE coordinates) accordingly.
#'
#' @param clustering_algorithm
#'   Character or `NULL`; name of the clustering algorithm to use for the
#'   prepared object (e.g., `"leiden"`, `"flowsom"`, `"louvain"`, `"phenograph"`).
#'   If `NULL` (default), uses the first detected algorithm.
#'   selected algorithm will be renamed to `"cluster"` for FCView compatibility.
#'
#' @param output_dir
#'   Character or `NULL`; path to directory where the .RData file will be saved.
#'   Must exist. If `NULL` (default), the object is not saved to disk.
#'
#' @param file_name
#'   Character or `NULL`; name of the .RData file (extension added automatically
#'   if not provided). If `NULL` (default), the object is not saved to disk.
#'   Both `output_dir` and `file_name` must be specified to enable saving.
#'
#' @param keep_fields
#'   Character vector; names of top‐level fields to retain in the prepared
#'   object. Default includes essential FCView fields: `"data"`, `"source"`,
#'   `"metadata"`, `"run_date"`, `"cluster"`, `"umap"`, `"tsne"`,
#'   `"cluster_heatmap"`, `"cluster_mapping"`, and `"metadata_sample"`.
#'   When `"cluster_mapping"` is included and the object contains a
#'   `{clustering_algorithm}_mapping` element (created by
#'   `FCSimple::fcs_annotate_clusters()`), that element is copied into the
#'   prepared object as `cluster_mapping`. If `"cluster_mapping"` is absent
#'   from `keep_fields`, no mapping element is transferred.
#'
#' @param scenith_compatible
#'   Logical; defaults to `FALSE`. When `TRUE`, the function expects a
#'   `scenith` slot in `fcs_join_obj` containing a data frame with one row
#'   per cell (in the same row order as `fcs_join_obj$data`), holding per-cell
#'   puromycin intensity values and associated metadata. This data frame is
#'   validated, optionally column-trimmed, and enriched with cluster assignments
#'   **before** any downsampling occurs, ensuring the full-resolution SCENITH
#'   data (all cells) is preserved in the prepared object for use in the FCView
#'   SCENITH tab.\cr\cr
#'   **Required columns in `obj$scenith`:**
#'   \itemize{
#'     \item `patient_ID` — sample identifier
#'     \item a puromycin intensity column matching the regex
#'       `^(puromycin|PURO|PUROMYCIN|Puromycin|Puro|puro)$`
#'     \item `inhibitor` — metabolic inhibitor condition label
#'   }
#'   **Optional columns:** any additional columns whose names also appear in
#'   `obj$metadata` (e.g., `timepoint`, `sample_ID`) are retained; all others
#'   are dropped. A `cluster` column is appended automatically from the
#'   selected clustering algorithm before downsampling.
#'
#' @param precalculate_scenith_centers
#'   Logical; defaults to `FALSE`. When `TRUE` and `scenith_compatible = TRUE`,
#'   the function pre-computes all four center estimators (mean, median, geometric mean,
#'   Hodges-Lehmann) for each inhibitor and stores them in a `scenith_centers` slot
#'   as a nested list: `scenith_centers[[inhibitor]][[estimator]]` = sample × cluster matrix.
#'   This front-loads computational cost at preparation time (typically 5-20 minutes per file,
#'   depending on cell count) but dramatically speeds up FCView app interactions on low-spec
#'   hardware, since the app skips center calculations and performs only higher-level
#'   operations (ratio, subtraction, etc.). If `FALSE` (default), the app computes centers
#'   on-the-fly as before. This parameter is ignored if `scenith_compatible = FALSE`.
#'
#' @param scenith_verbose
#'   Logical; defaults to `TRUE`. When `TRUE`, the function emits progress
#'   messages while precalculating SCENITH centers, including inhibitor, sample,
#'   and cluster-level progress.
#'
#' @param scenith_sample_col
#'   Character; the sample unit column name for SCENITH precalculation.
#'   Defaults to "patient_ID".
#'
#' @param precalculate_mfi
#'   Logical; defaults to `TRUE`. When `TRUE`, the function pre-computes
#'   per-sample, per-cluster marker intensity center matrices for each marker
#'   in the expression data. Computations run on the full (undownsampled)
#'   dataset, mirroring how `precalculate_scenith_centers` handles SCENITH
#'   data. The result is stored as `obj$mfi`, a nested list with structure
#'   `mfi[[statistic]][[marker]]`, where each element is a samples × clusters
#'   numeric matrix. Supported statistics are controlled by `mfi_stats`.
#'
#' @param mfi_markers
#'   Character vector or `NULL`; which expression matrix columns (markers) to
#'   precompute MFI matrices for. If `NULL` (default), all columns in
#'   `fcs_join_obj$data` are used. Invalid marker names trigger an error.
#'
#' @param mfi_stats
#'   Character vector; which center statistics to precompute. Must be a subset
#'   of `c("mean", "median")`. Defaults to both. Geometric mean is intentionally
#'   excluded — it cannot be reliably re‑aggregated from per‑cluster values when
#'   annotations combine multiple clusters into celltypes.
#'
#' @details
#'   The function performs the following steps:
#'   1. Validates required fields (`data`, `source`, `metadata`) are present.
#'   2. Detects all clustering algorithms in the object.
#'   3. Selects the specified or first detected algorithm.
#'   4. Optionally downsamples cells (with fixed seed 123 for reproducibility),
#'      updating all cell‐level fields: data matrix, source, run_date,
#'      cluster assignments, and UMAP/tSNE coordinates.
#'   5. Renames the selected clustering algorithm to `"cluster"` and its heatmap
#'      to `"cluster_heatmap"` (required for FCView app compatibility).
#'      The existing `frequency` and `counts` tables are preserved as-is from
#'      the original full-dataset calculation; they are **never** recalculated.
#'      An error is raised if either table is absent — ensure
#'      `FCSimple::fcs_calculate_abundance()` has been run beforehand.
#'   6. Removes all other clustering algorithms and their heatmaps.
#'   7. Filters the object to retain only fields specified in `keep_fields`.
#'      If `"cluster_mapping"` is in `keep_fields`, the function additionally
#'      resolves the corresponding `{selected_algo}_mapping` element (if present)
#'      and includes it as `cluster_mapping` in the prepared object.
#'   8. Cleans up cluster and heatmap structures to keep only essential elements:
#'      - For `cluster`: keeps `clusters`, `settings`, `frequency`, `fraction`, `counts`
#'      - For `cluster_heatmap`: keeps `heatmap_tile_data`, `population_size`, `rep_used`
#'   9. Ensures data types match FCView expectations (data as matrix, coordinates
#'      as data frames).
#'   10. Adds metadata tracking: `fcview_prepared`, `fcview_prep_time`,
#'       `fcview_clustering_algorithm`, and optionally `fcview_downsampled`.
#'   11. If `output_dir` and `file_name` are provided, saves the object as
#'       `fcs_data` in an .RData file and reports file size.
#'
#' @return
#'   The prepared `fcs_join_obj`, with non‐essential fields removed, the
#'   selected clustering algorithm renamed to `"cluster"`, and metadata
#'   tracking added. The object is returned invisibly if saved to disk,
#'   otherwise returned explicitly.
#'
#' @examples
#' \dontrun{
#'   # Basic preparation without saving (use packaged example files)
#'   files <- FCSimple::fcs_example_files()
#'   joined <- FCSimple::fcs_join(files)
#'   clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#'   reduced <- FCSimple::fcs_reduce_dimensions(clustered, algorithm = "umap")
#'
#'   prepared <- FCSimple::fcs_prepare_fcview_object(reduced)
#'
#'   # Prepare with downsampling to 50,000 cells
#'   prepared_ds <- FCSimple::fcs_prepare_fcview_object(
#'     reduced,
#'     downsample_size = 50000
#'   )
#'
#'   # Prepare and save to file
#'   prepared_saved <- FCSimple::fcs_prepare_fcview_object(
#'     reduced,
#'     downsample_size = 100000,
#'     clustering_algorithm = "leiden",
#'     output_dir = "output/fcview",
#'     file_name = "my_analysis_fcview"
#'   )
#'
#'   # Prepare with specific algorithm selection
#'   prepared_flowsom <- FCSimple::fcs_prepare_fcview_object(
#'     reduced,
#'     clustering_algorithm = "flowsom",
#'     output_dir = getwd(),
#'     file_name = "flowsom_analysis.RData"
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, FCSimple::fcs_cluster, FCSimple::fcs_reduce_dimensions,
#'   FCSimple::fcs_cluster_heatmap
#'
#' @export
fcs_prepare_fcview_object <- function(fcs_join_obj,
                                      downsample_size = NULL,
                                      clustering_algorithm = NULL,
                                      output_dir = NULL,
                                      file_name = NULL,
                                      scenith_compatible = FALSE,
                                      precalculate_scenith_centers = FALSE,
                                      scenith_sample_col = "patient_ID",
                                      scenith_verbose = TRUE,
                                      precalculate_mfi = TRUE,
                                      mfi_markers = NULL,
                                      mfi_stats = c("mean", "median"),
                                      keep_fields = c("data", "source", "metadata", "run_date",
                                                      "cluster", "umap", "tsne",
                                                      "cluster_heatmap", "cluster_mapping",
                                                      "metadata_sample")) {

  if (!is.list(fcs_join_obj)) {
    stop("fcs_join_obj must be a list")
  }

  if (!'metadata' %in% names(fcs_join_obj)) {
    meta <- data.frame(patient_ID = fcs_join_obj$source, run_date = fcs_join_obj$run_date)
    meta <- meta[!duplicated(meta$patient_ID),]
    fcs_join_obj$metadata <- meta
  }

  required_fields <- c("data", "source", "metadata")
  missing_required <- setdiff(required_fields, names(fcs_join_obj))
  if (length(missing_required) > 0) {
    stop("Missing required fields: ", paste(missing_required, collapse = ", "))
  }

  all_cluster_algos <- c("leiden", "flowsom", "louvain", "phenograph", "cluster")
  present_algos <- c()
  for (algo in all_cluster_algos) {
    if (algo %in% names(fcs_join_obj)) {
      if (is.list(fcs_join_obj[[algo]]) && "clusters" %in% names(fcs_join_obj[[algo]])) {
        present_algos <- c(present_algos, algo)
      }
    }
  }

  if (length(present_algos) == 0) {
    stop("No clustering algorithm found in fcs_join_obj. Expected one of: ",
         paste(all_cluster_algos, collapse = ", "))
  }

  if (!is.null(clustering_algorithm)) {
    if (!clustering_algorithm %in% present_algos) {
      stop("Specified clustering_algorithm '", clustering_algorithm,
           "' not found in object. Available: ", paste(present_algos, collapse = ", "))
    }
    selected_algo <- clustering_algorithm
  } else {
    selected_algo <- present_algos[1]
    if (length(present_algos) > 1) {
      message(paste0("Multiple clustering algorithms detected: ",
                     paste(present_algos, collapse = ", ")))
      message(paste0("Using '", selected_algo, "' (first detected). ",
                     "Specify clustering_algorithm to choose a different one."))
    }
  }

  n_cells <- nrow(fcs_join_obj$data)
  n_cells_original <- n_cells

  # ---- SCENITH: validate and enrich SCENITH data BEFORE downsampling ----
  if (scenith_compatible) {
    if (!"scenith" %in% names(fcs_join_obj)) {
      stop("scenith_compatible = TRUE but fcs_join_obj$scenith was not found. ",
           "Ensure the object contains a 'scenith' data.frame with one row per cell.")
    }
    puro_df <- fcs_join_obj$scenith
    if (!is.data.frame(puro_df)) {
      stop("fcs_join_obj$scenith must be a data.frame.")
    }
    if (nrow(puro_df) != n_cells) {
      stop("nrow(fcs_join_obj$scenith) (", nrow(puro_df), ") must equal ",
           "nrow(fcs_join_obj$data) (", n_cells, "). ",
           "Ensure the scenith data.frame has one row per cell in the same order as the data matrix.")
    }

    # Validate `scenith_sample_col`
    if (!scenith_sample_col %in% colnames(puro_df)) {
      stop("The specified scenith_sample_col ('", scenith_sample_col, "') does not exist in the scenith data frame.")
    }

    # Required column: puromycin intensity (flexible naming)
    puro_col_pattern <- "^(puromycin|PURO|PUROMYCIN|Puromycin|Puro|puro)$"
    puro_col <- grep(puro_col_pattern, colnames(puro_df), value = TRUE)
    if (length(puro_col) == 0) {
      stop("fcs_join_obj$scenith must contain a puromycin intensity column matching ",
           "the regex ^(puromycin|PURO|PUROMYCIN|Puromycin|Puro|puro)$")
    }
    if (length(puro_col) > 1) {
      warning("Multiple puromycin intensity columns found; using '", puro_col[1], "'.")
      puro_col <- puro_col[1]
    }

    # Required column: inhibitor
    if (!"inhibitor" %in% colnames(puro_df)) {
      stop("fcs_join_obj$scenith must contain an 'inhibitor' column.")
    }

    # Determine columns to keep: required only.
    final_puro_cols <- c(scenith_sample_col, puro_col, "inhibitor")

    # Verify that every sample unit in scenith exists in metadata
    if (is.data.frame(fcs_join_obj$metadata) && scenith_sample_col %in% colnames(fcs_join_obj$metadata)) {
      puro_ids <- unique(puro_df[[scenith_sample_col]])
      meta_ids <- unique(fcs_join_obj$metadata[[scenith_sample_col]])
      missing_ids <- setdiff(puro_ids, meta_ids)
      if (length(missing_ids) > 0) {
        stop("The following ", scenith_sample_col, " value(s) appear in fcs_join_obj$scenith but not in ",
             "fcs_join_obj$metadata: ", paste(missing_ids, collapse = ", "),
             ". Ensure both data sources share the same ", scenith_sample_col, " values.")
      }
    } else {
      warning("fcs_join_obj$metadata does not contain the specified scenith_sample_col ('", scenith_sample_col, "'); ",
              "cannot verify scenith / metadata alignment.")
    }

    # Add cluster assignments from the full (pre-downsample) dataset
    cluster_vec <- fcs_join_obj[[selected_algo]]$clusters
    if (length(cluster_vec) != n_cells) {
      stop("Length of cluster assignments (", length(cluster_vec), ") does not match ",
           "nrow(fcs_join_obj$data) (", n_cells, ").")
    }
    puro_df$cluster <- cluster_vec
    final_puro_cols <- c(final_puro_cols, "cluster")

    fcs_join_obj$scenith <- puro_df[, final_puro_cols, drop = FALSE]

    # NA check — no NAs are permitted in any column
    na_counts <- colSums(is.na(fcs_join_obj$scenith))
    cols_with_na <- names(na_counts[na_counts > 0])
    if (length(cols_with_na) > 0) {
      stop("fcs_join_obj$scenith contains NA values in the following column(s): ",
           paste(sprintf("'%s' (%d NA%s)", cols_with_na, na_counts[cols_with_na],
                         ifelse(na_counts[cols_with_na] == 1, "", "s")), collapse = ", "),
           ". Remove or impute NA values before preparing the FCView object.")
    }

    message("SCENITH: data validated and cluster assignments added (pre-downsample).")
    message(paste0("  SCENITH rows (all cells, not downsampled): ", nrow(fcs_join_obj$scenith)))
    message(paste0("  Puromycin intensity column: '", puro_col, "'"))
    message("  Columns: ", scenith_sample_col, ", ", puro_col, ", inhibitor, cluster")

    # ---- Precalculate SCENITH centers if requested ----
    if (precalculate_scenith_centers) {
      message("[fcs_prepare_fcview_object] Precalculating SCENITH centers...")
      if (isTRUE(scenith_verbose)) {
        message("[fcs_prepare_fcview_object] SCENITH precalc sample column: '",
                scenith_sample_col, "'")
      }

      # Helper function to compute all four center estimators for an inhibitor
      compute_inhibitor_centers <- function(scenith_df, puro_col_name) {
        # scenith_df: filtered data.frame for one inhibitor, with columns:
        #   scenith_sample_col, {puro_col_name}, cluster
        # Returns: list(mean, median, geomean, hl), each a sample x cluster matrix

        samples <- unique(scenith_df[[scenith_sample_col]])
        clusters <- unique(scenith_df$cluster)

        if (isTRUE(scenith_verbose)) {
          message("[fcs_prepare_fcview_object]   ", length(samples), " sample unit(s), ",
                  length(clusters), " cluster(s)")
        }

        # Initialize result matrices
        centers <- list(
          mean    = matrix(NA, nrow = length(samples), ncol = length(clusters),
                           dimnames = list(samples, as.character(clusters))),
          median  = matrix(NA, nrow = length(samples), ncol = length(clusters),
                           dimnames = list(samples, as.character(clusters))),
          geomean = matrix(NA, nrow = length(samples), ncol = length(clusters),
                           dimnames = list(samples, as.character(clusters))),
          hl      = matrix(NA, nrow = length(samples), ncol = length(clusters),
                           dimnames = list(samples, as.character(clusters)))
        )

        for (sample in samples) {
          sample_df <- scenith_df[scenith_df[[scenith_sample_col]] == sample, ]

          if (isTRUE(scenith_verbose)) {
            message("[fcs_prepare_fcview_object]   sample '", sample, "' (", nrow(sample_df),
                    " row(s))")
          }

          for (cluster_id in clusters) {
            entity_vals <- sample_df[[puro_col_name]][sample_df$cluster == cluster_id]

            if (length(entity_vals) > 0) {
              if (isTRUE(scenith_verbose)) {
                message("[fcs_prepare_fcview_object]     cluster '", cluster_id, "': ",
                        length(entity_vals), " cell(s)")
              }

              # Mean
              centers$mean[sample, as.character(cluster_id)] <- mean(entity_vals, na.rm = TRUE)

              # Median
              centers$median[sample, as.character(cluster_id)] <- median(entity_vals, na.rm = TRUE)

              # Geometric mean
              valid_vals <- entity_vals[entity_vals > 0]
              if (length(valid_vals) > 0) {
                centers$geomean[sample, as.character(cluster_id)] <- exp(mean(log(valid_vals)))
              }

              # Hodges-Lehmann (approximate via Wilcoxon)
              if (length(entity_vals) >= 2) {
                tryCatch({
                  # Subsample to 5000 cells for speed (same as app.R logic)
                  HL_N_MAX <- 5000L
                  entity_vals_hl <- if (length(entity_vals) > HL_N_MAX) {
                    set.seed(123); sample(entity_vals, HL_N_MAX)
                  } else {
                    entity_vals
                  }
                  centers$hl[sample, as.character(cluster_id)] <-
                    as.numeric(wilcox.test(entity_vals_hl, conf.int = TRUE, exact = FALSE)$estimate)
                }, error = function(e) {
                  # If Wilcoxon fails, leave as NA
                  NA
                })
              }
            }
          }
        }

        centers
      }

      # Compute centers for all inhibitors
      scenith_centers <- list()
      inhibitors <- unique(fcs_join_obj$scenith$inhibitor)

      if (isTRUE(scenith_verbose)) {
        message("[fcs_prepare_fcview_object] Inhibitors to process: ",
                paste(inhibitors, collapse = ", "))
      }

      for (inh_idx in seq_along(inhibitors)) {
        inh <- inhibitors[[inh_idx]]
        if (isTRUE(scenith_verbose)) {
          message("[fcs_prepare_fcview_object] Processing inhibitor ", inh_idx, "/",
                  length(inhibitors), ": '", inh, "'")
        }
        inh_df <- fcs_join_obj$scenith[fcs_join_obj$scenith$inhibitor == inh, ]
        scenith_centers[[inh]] <- compute_inhibitor_centers(inh_df, puro_col)
        if (isTRUE(scenith_verbose)) {
          message("[fcs_prepare_fcview_object] Finished inhibitor '", inh, "'")
        }
      }

      # Store in the object (will be included if keep_fields contains "scenith_centers")
      fcs_join_obj$scenith_centers <- scenith_centers
      fcs_join_obj$scenith_centers_sample_col <- scenith_sample_col
      keep_fields <- union(keep_fields, "scenith_centers")
      keep_fields <- union(keep_fields, "scenith_centers_sample_col")

      message("[fcs_prepare_fcview_object] SCENITH centers precalculation complete.")
    }
  }
  # ---- end SCENITH block ----

  # ---- MFI: precalculate per-marker sample × cluster center matrices BEFORE downsampling ----
  if (precalculate_mfi) {
    mfi_stats <- match.arg(mfi_stats, choices = c("mean", "median"), several.ok = TRUE)
    if (length(mfi_stats) == 0) {
      stop("mfi_stats must contain at least one of 'mean' or 'median'")
    }

    # Determine which markers to compute (from data slot)
    if (is.null(mfi_markers)) {
      mfi_markers <- colnames(fcs_join_obj$data)
    } else {
      missing_markers <- setdiff(mfi_markers, colnames(fcs_join_obj$data))
      if (length(missing_markers) > 0) {
        stop("The following mfi_markers are not columns in fcs_join_obj$data: ",
             paste(missing_markers, collapse = ", "))
      }
    }

    # Get source and cluster vectors (full, pre-downsample)
    source_vec <- fcs_join_obj$source
    cluster_vec <- fcs_join_obj[[selected_algo]]$clusters
    samples <- unique(source_vec)
    clusters <- sort(unique(cluster_vec))
    n_samples <- length(samples)
    n_clusters <- length(clusters)

    # ---- Helper: compute MFI for a given expression matrix and a set of markers ----
    compute_mfi_layer <- function(expr_mat, markers, source_vec, cluster_vec,
                                  mfi_stats, samples, clusters, n_samples, n_clusters,
                                  layer_label) {
      message(sprintf("[fcs_prepare_fcview_object] Computing MFI from %s data...", layer_label))
      message("  Markers: ", length(markers))
      message("  Statistics: ", paste(mfi_stats, collapse = ", "))

      mfi_layer <- list()
      for (stat in mfi_stats) {
        mfi_layer[[stat]] <- list()
      }

      for (marker in markers) {
        marker_vals <- expr_mat[, marker]

        df <- data.frame(
          sample  = source_vec,
          cluster = as.character(cluster_vec),
          value   = as.numeric(marker_vals),
          stringsAsFactors = FALSE
        )

        # ---- Mean ----
        if ("mean" %in% mfi_stats) {
          sum_mat <- xtabs(value ~ sample + cluster, data = df)
          cnt_mat <- xtabs(~ sample + cluster, data = df)
          cnt_mat[cnt_mat == 0] <- NA
          mean_mat <- sum_mat / cnt_mat
          storage.mode(mean_mat) <- "double"

          mean_mat_full <- matrix(NA_real_, nrow = n_samples, ncol = n_clusters,
                                  dimnames = list(samples, as.character(clusters)))
          for (s in rownames(mean_mat)) {
            for (c in colnames(mean_mat)) {
              if (!is.na(mean_mat[s, c])) {
                mean_mat_full[s, c] <- mean_mat[s, c]
              }
            }
          }
          mfi_layer[["mean"]][[marker]] <- mean_mat_full
        }

        # ---- Median ----
        if ("median" %in% mfi_stats) {
          med_agg <- aggregate(value ~ sample + cluster, data = df, FUN = median, na.rm = TRUE)
          median_mat_full <- matrix(NA_real_, nrow = n_samples, ncol = n_clusters,
                                    dimnames = list(samples, as.character(clusters)))
          for (r in seq_len(nrow(med_agg))) {
            s <- as.character(med_agg$sample[r])
            c <- as.character(med_agg$cluster[r])
            median_mat_full[s, c] <- med_agg$value[r]
          }
          mfi_layer[["median"]][[marker]] <- median_mat_full
        }
      }

      message(sprintf("[fcs_prepare_fcview_object] %s data MFI complete: %d marker(s) × %d statistic(s)",
                      layer_label, length(markers), length(mfi_stats)))
      mfi_layer
    }

    # ---- Compute MFI from transformed (data) slot ----
    mfi_list <- list()
    message("[fcs_prepare_fcview_object] Precalculating MFI matrices...")
    message("  Computation on full (undownsampled) dataset: ", n_cells, " cells")

    mfi_list[["data"]] <- compute_mfi_layer(
      expr_mat    = fcs_join_obj$data,
      markers     = mfi_markers,
      source_vec  = source_vec,
      cluster_vec = cluster_vec,
      mfi_stats   = mfi_stats,
      samples     = samples,
      clusters    = clusters,
      n_samples   = n_samples,
      n_clusters  = n_clusters,
      layer_label = "transformed"
    )

    # ---- Conditionally compute MFI from raw slot ----
    has_raw <- !is.null(fcs_join_obj$raw) &&
               !identical(fcs_join_obj$raw, NA) &&
               is.matrix(fcs_join_obj$raw)

    if (has_raw) {
      # Determine if raw is actually different from data (transformations were applied)
      raw_differs <- !identical(fcs_join_obj$raw, fcs_join_obj$data)

      if (raw_differs) {
        # Use intersection of markers present in both raw and data
        raw_markers <- intersect(mfi_markers, colnames(fcs_join_obj$raw))
        if (length(raw_markers) > 0) {
          message("[fcs_prepare_fcview_object] Raw data slot detected and differs from data. ",
                  "Computing MFI for ", length(raw_markers), " shared marker(s).")
          mfi_list[["raw"]] <- compute_mfi_layer(
            expr_mat    = fcs_join_obj$raw,
            markers     = raw_markers,
            source_vec  = source_vec,
            cluster_vec = cluster_vec,
            mfi_stats   = mfi_stats,
            samples     = samples,
            clusters    = clusters,
            n_samples   = n_samples,
            n_clusters  = n_clusters,
            layer_label = "untransformed (raw)"
          )
        } else {
          message("[fcs_prepare_fcview_object] Raw data slot present but no shared markers with data. ",
                  "Skipping raw MFI computation.")
        }
      }
    }

    # Store in the object
    fcs_join_obj$mfi <- mfi_list
    keep_fields <- union(keep_fields, "mfi")

    message("[fcs_prepare_fcview_object] MFI precalculation complete.")
    message("  Sources: ", paste(names(mfi_list), collapse = ", "))
    for (src_name in names(mfi_list)) {
      n_markers <- length(mfi_list[[src_name]][["mean"]])
      message(sprintf("  [%s] %d marker(s) × %d statistic(s), %d samples × %d clusters",
                      src_name, n_markers, length(mfi_stats), n_samples, n_clusters))
    }
  }
  # ---- end MFI block ----

  if (length(fcs_join_obj$source) != n_cells) {
    stop("source length must equal nrow(data)")
  }

  # run_date must be either absent or cell-length — a scalar survives upload
  # but is rejected by FCView's validateInput.
  if ("run_date" %in% names(fcs_join_obj) &&
      !is.null(fcs_join_obj$run_date) &&
      length(fcs_join_obj$run_date) != n_cells) {
    stop("run_date length (", length(fcs_join_obj$run_date), ") must equal nrow(data) (", n_cells, "). ",
         "Ensure run_date is a per-cell vector, not a scalar or per-sample summary.")
  }


  if (!is.null(downsample_size)) {
    if (!is.numeric(downsample_size) || downsample_size < 1 || downsample_size %% 1 != 0) {
      stop("downsample_size must be a positive integer")
    }

    if (downsample_size >= n_cells) {
      warning("downsample_size is >= number of cells. No downsampling performed.")
    } else {
      message(paste0("Downsampling from ", n_cells, " to ", downsample_size, " cells..."))
      set.seed(123)
      keep_idx <- sample(1:n_cells, size = downsample_size, replace = FALSE)

      fcs_join_obj$data <- fcs_join_obj$data[keep_idx, , drop = FALSE]
      fcs_join_obj$source <- fcs_join_obj$source[keep_idx]

      if ("run_date" %in% names(fcs_join_obj) && length(fcs_join_obj$run_date) == n_cells) {
        fcs_join_obj$run_date <- fcs_join_obj$run_date[keep_idx]
      }

      if (selected_algo %in% names(fcs_join_obj)) {
        if (is.list(fcs_join_obj[[selected_algo]]) && "clusters" %in% names(fcs_join_obj[[selected_algo]])) {
          if (length(fcs_join_obj[[selected_algo]]$clusters) == n_cells) {
            fcs_join_obj[[selected_algo]]$clusters <- fcs_join_obj[[selected_algo]]$clusters[keep_idx]
          }
        }
      }

      if ("umap" %in% names(fcs_join_obj)) {
        if (is.list(fcs_join_obj$umap) && "coordinates" %in% names(fcs_join_obj$umap)) {
          if (is.data.frame(fcs_join_obj$umap$coordinates) || is.matrix(fcs_join_obj$umap$coordinates)) {
            if (nrow(fcs_join_obj$umap$coordinates) == n_cells) {
              fcs_join_obj$umap$coordinates <- fcs_join_obj$umap$coordinates[keep_idx, , drop = FALSE]
            }
          }
        }
      }

      if ("tsne" %in% names(fcs_join_obj)) {
        if (is.list(fcs_join_obj$tsne) && "coordinates" %in% names(fcs_join_obj$tsne)) {
          if (is.data.frame(fcs_join_obj$tsne$coordinates) || is.matrix(fcs_join_obj$tsne$coordinates)) {
            if (nrow(fcs_join_obj$tsne$coordinates) == n_cells) {
              fcs_join_obj$tsne$coordinates <- fcs_join_obj$tsne$coordinates[keep_idx, , drop = FALSE]
            }
          }
        }
      }

      n_cells <- downsample_size
    }
  }

  # Rename selected algorithm to cluster (no recalculation — preserve original

  # full-dataset frequency and counts tables as-is)
  if (selected_algo != "cluster") {
    fcs_join_obj$cluster <- fcs_join_obj[[selected_algo]]
    selected_heatmap <- paste0(selected_algo, "_heatmap")
    if (selected_heatmap %in% names(fcs_join_obj)) {
      fcs_join_obj$cluster_heatmap <- fcs_join_obj[[selected_heatmap]]
    }
  }

  if (!"frequency" %in% names(fcs_join_obj$cluster)) {
    stop("cluster$frequency is missing. Run fcs_calculate_abundance() before preparing the FCView object.")
  }
  if (!"counts" %in% names(fcs_join_obj$cluster)) {
    stop("cluster$counts is missing. Run fcs_calculate_abundance() before preparing the FCView object.")
  }

  algos_to_remove <- setdiff(all_cluster_algos, "cluster")
  for (algo in algos_to_remove) {
    if (algo %in% names(fcs_join_obj)) {
      fcs_join_obj[[algo]] <- NULL
    }
    heatmap_name <- paste0(algo, "_heatmap")
    if (heatmap_name %in% names(fcs_join_obj)) {
      fcs_join_obj[[heatmap_name]] <- NULL
    }
  }

  # Ensure scenith is included when scenith_compatible = TRUE
  if (scenith_compatible) {
    keep_fields <- union(keep_fields, "scenith")
  }

  prepared_obj <- list()
  for (field in names(fcs_join_obj)) {
    if (field %in% keep_fields) {
      prepared_obj[[field]] <- fcs_join_obj[[field]]
    }
  }

  # Resolve {algo}_mapping → cluster_mapping when requested
  if ("cluster_mapping" %in% keep_fields) {
    mapping_name <- paste0(selected_algo, "_mapping")
    if (mapping_name %in% names(fcs_join_obj)) {
      prepared_obj$cluster_mapping <- fcs_join_obj[[mapping_name]]
    }
    # If no {algo}_mapping exists (user never called fcs_annotate_clusters),
    # cluster_mapping is simply absent from the prepared object — FCView
    # will behave as if no pre-annotation was provided.
  }

  if (!is.null(prepared_obj$data) && !is.matrix(prepared_obj$data)) {
    prepared_obj$data <- as.matrix(prepared_obj$data)
  }

  if ("umap" %in% names(prepared_obj) && !is.null(prepared_obj$umap$coordinates)) {
    if (is.matrix(prepared_obj$umap$coordinates)) {
      prepared_obj$umap$coordinates <- as.data.frame(prepared_obj$umap$coordinates)
    }
  }

  if ("tsne" %in% names(prepared_obj) && !is.null(prepared_obj$tsne$coordinates)) {
    if (is.matrix(prepared_obj$tsne$coordinates)) {
      prepared_obj$tsne$coordinates <- as.data.frame(prepared_obj$tsne$coordinates)
    }
  }

  if (is.data.frame(prepared_obj$metadata)) {
    if (!"patient_ID" %in% colnames(prepared_obj$metadata)) {
      warning("metadata does not contain 'patient_ID' column. Consider adding it for app.R compatibility.")
    }
  }

  if ("cluster" %in% names(prepared_obj)) {
    if (is.list(prepared_obj$cluster)) {
      essential_cluster_fields <- c("clusters", "settings", "frequency", "fraction", "counts")
      cluster_fields_present <- intersect(essential_cluster_fields, names(prepared_obj$cluster))
      new_cluster <- list()
      for (cf in cluster_fields_present) {
        new_cluster[[cf]] <- prepared_obj$cluster[[cf]]
      }
      prepared_obj$cluster <- new_cluster
    }
  }

  if ("cluster_heatmap" %in% names(prepared_obj) && is.list(prepared_obj$cluster_heatmap)) {
    essential_heatmap_fields <- c("heatmap_tile_data", "population_size", "rep_used")
    heatmap_fields_present <- intersect(essential_heatmap_fields, names(prepared_obj$cluster_heatmap))
    new_heatmap <- list()
    for (hf in heatmap_fields_present) {
      new_heatmap[[hf]] <- prepared_obj$cluster_heatmap[[hf]]
    }
    prepared_obj$cluster_heatmap <- new_heatmap
  }

  prepared_obj$fcview_prepared <- TRUE
  prepared_obj$fcview_prep_time <- Sys.time()
  prepared_obj$fcview_clustering_algorithm <- selected_algo
  if (!is.null(downsample_size) && downsample_size < n_cells_original) {
    prepared_obj$fcview_downsampled <- TRUE
    prepared_obj$fcview_original_ncells <- n_cells_original
    prepared_obj$fcview_final_ncells <- n_cells
  }

  message("FCView object prepared successfully.")
  message(paste0("  Cells: ", n_cells))
  message(paste0("  Features: ", ncol(prepared_obj$data)))
  message(paste0("  Clustering algorithm: ", selected_algo, " (renamed to 'cluster')"))

  if ("umap" %in% names(prepared_obj)) message("  UMAP coordinates: present")
  if ("tsne" %in% names(prepared_obj)) message("  tSNE coordinates: present")
  if ("cluster_heatmap"  %in% names(prepared_obj)) message("  Cluster heatmap: present")
  if ("cluster_mapping"  %in% names(prepared_obj)) {
    n_mapped <- sum(!is.na(prepared_obj$cluster_mapping$celltype))
    n_total_m <- nrow(prepared_obj$cluster_mapping)
    message(sprintf("  Cluster mapping: present (%d/%d clusters annotated, from '%s')",
                    n_mapped, n_total_m, paste0(selected_algo, "_mapping")))
  }
  if (scenith_compatible && "scenith" %in% names(prepared_obj)) {
    message(paste0("  SCENITH data: present (", nrow(prepared_obj$scenith),
                   " rows \u00d7 ", ncol(prepared_obj$scenith), " columns, undownsampled)"))
    if ("scenith_centers" %in% names(prepared_obj)) {
      n_inh <- length(prepared_obj$scenith_centers)
      message(paste0("  SCENITH precalc centers: present (", n_inh, " inhibitor(s), 4 estimators each)"))
      if ("scenith_centers_sample_col" %in% names(prepared_obj)) {
        message(paste0("  SCENITH precalc sample column: '", prepared_obj$scenith_centers_sample_col, "'"))
      }
    }
  }
  if ("mfi" %in% names(prepared_obj)) {
    sources <- names(prepared_obj$mfi)
    for (src_name in sources) {
      n_markers <- length(prepared_obj$mfi[[src_name]][["mean"]])
      n_stats <- length(prepared_obj$mfi[[src_name]])
      stat_names <- paste(names(prepared_obj$mfi[[src_name]]), collapse = ", ")
      source_label <- if (src_name == "data") "Transformed" else if (src_name == "raw") "Untransformed" else src_name
      message(paste0("  MFI [", source_label, "]: present (", n_markers, " marker(s), ",
                     n_stats, " statistic(s): ", stat_names, ")"))
    }
  }

  # Save to file if output_dir and file_name are provided
  if (!is.null(output_dir) && !is.null(file_name)) {
    # Validate output directory
    if (!dir.exists(output_dir)) {
      stop("Output directory does not exist: ", output_dir)
    }

    # Ensure .RData extension
    if (!grepl("\\.RData$", file_name, ignore.case = TRUE)) {
      file_name <- paste0(file_name, ".RData")
    }

    # Construct full path
    full_path <- file.path(output_dir, file_name)

    # Save the object
    message("\nSaving FCView object...")
    fcs_data <- prepared_obj
    save(fcs_data, file = full_path)

    # Calculate file size
    file_size_mb <- round(file.info(full_path)$size / (1024^2), 2)

    message(paste0("  File saved: ", full_path))
    message(paste0("  File size: ", file_size_mb, " MB"))
  } else if (!is.null(output_dir) || !is.null(file_name)) {
    warning("Both output_dir and file_name must be specified to save the object. Skipping save.")
  }

  return(prepared_obj)
}
