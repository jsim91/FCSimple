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
#'
#' @param downsample_size
#'   Integer or `NULL`; maximum number of cells to retain. If `NULL` (default),
#'   no downsampling is performed. If specified, randomly samples cells and
#'   updates all associated fields (data, source, run_date, cluster assignments,
#'   UMAP/tSNE coordinates) accordingly.
#'
#' @param clustering_algorithm
#'   Character or `NULL`; name of the clustering algorithm to use for the
#'   prepared object (e.g., `"leiden"`, `"flowsom"`, `"louvain"`, `"phenograph"`,
#'   `"git"`). If `NULL` (default), uses the first detected algorithm. The
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
#'   6. Removes all other clustering algorithms and their heatmaps.
#'   7. Filters the object to retain only fields specified in `keep_fields`.
#'   8. Cleans up cluster and heatmap structures to keep only essential elements:
#'      - For `cluster`: keeps `clusters`, `settings`, `abundance`, `counts`
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
#'   # Basic preparation without saving
#'   files <- list.files("data/fcs", pattern = "\\.fcs$", full.names = TRUE)
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

  all_cluster_algos <- c("leiden", "flowsom", "louvain", "phenograph", "git", "cluster")
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

  if (length(fcs_join_obj$source) != n_cells) {
    stop("source length must equal nrow(data)")
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

  # Rename selected algorithm to cluster
  if (selected_algo != "cluster") {
    fcs_join_obj$cluster <- fcs_join_obj[[selected_algo]]
    if('abundance' %in% names(fcs_join_obj$cluster)) {
      fcs_join_obj <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_join_obj,
                                                        report_algorithm = 'cluster',
                                                        report_as = 'frequency')
    }
    if('counts' %in% names(fcs_join_obj$cluster)) {
      fcs_join_obj <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_join_obj,
                                                        report_algorithm = 'cluster',
                                                        report_as = 'count')
    }
    selected_heatmap <- paste0(selected_algo, "_heatmap")
    if (selected_heatmap %in% names(fcs_join_obj)) {
      fcs_join_obj$cluster_heatmap <- fcs_join_obj[[selected_heatmap]]
    }
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

  prepared_obj <- list()
  for (field in names(fcs_join_obj)) {
    if (field %in% keep_fields) {
      prepared_obj[[field]] <- fcs_join_obj[[field]]
    }
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
      essential_cluster_fields <- c("clusters", "settings", "abundance", "counts")
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
  if (!is.null(downsample_size) && downsample_size < nrow(fcs_join_obj$data)) {
    prepared_obj$fcview_downsampled <- TRUE
    prepared_obj$fcview_original_ncells <- nrow(fcs_join_obj$data)
    prepared_obj$fcview_final_ncells <- n_cells
  }

  message("FCView object prepared successfully.")
  message(paste0("  Cells: ", n_cells))
  message(paste0("  Features: ", ncol(prepared_obj$data)))
  message(paste0("  Clustering algorithm: ", selected_algo, " (renamed to 'cluster')"))

  if ("umap" %in% names(prepared_obj)) message("  UMAP coordinates: present")
  if ("tsne" %in% names(prepared_obj)) message("  tSNE coordinates: present")
  if ("cluster_heatmap" %in% names(prepared_obj)) message("  Cluster heatmap: present")

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
