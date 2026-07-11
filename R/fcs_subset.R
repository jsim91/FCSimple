#' @title Subset FCSimple Analysis Object by Sample or Cluster
#'
#' @description
#'   Creates a new FCSimple-style object containing only the events that match
#'   specified sample names or cluster IDs. Can optionally retain per-cell
#'   slots such as clustering assignments, dimension reduction coordinates,
#'   and custom user vectors via the `keep_slots` argument. Records the
#'   subsetting criteria and original indices in a `subset_on` element for
#'   provenance.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), optionally processed by
#'   FCSimple::fcs_cluster(), FCSimple::fcs_reduce_dimensions(), or other
#'   functions. Should include at least `data` (numeric matrix) and
#'   `source` (vector of sample identifiers).
#'
#' @param subset_by
#'   Character; one of `"source"` or `"cluster"`. Determines whether to subset
#'   by sample names (`source`) or by cluster membership (`cluster`).
#'   Default `c("cluster","source")`.
#'
#' @param subset_cluster_algorithm
#'   Character or `NA`; when `subset_by = "cluster"`, the name of the clustering
#'   element in `fcs_join_obj` (e.g. `"leiden"`, `"flowsom"`). Ignored if
#'   `subset_by = "source"`. Default `NA`.
#'
#' @param subset_values
#'   Vector; the values to retain. If `subset_by = "source"`, a character
#'   vector of sample names. If `subset_by = "cluster"`, a numeric or character
#'   vector of cluster IDs.
#'
#' @param keep_slots
#'   Character vector or `NULL` (default). Names of additional slots in
#'   `fcs_join_obj` to subset and retain in the output. Supported slot types:
#'   \itemize{
#'     \item Clustering results: `"leiden"`, `"louvain"`, `"flowsom"`,
#'           `"phenograph"` — subsets the `$clusters` vector.
#'     \item Dimension reduction results: `"umap"`, `"tsne"` — subsets the
#'           `$coordinates` matrix.
#'     \item Simple per-cell vectors: any top-level element that is an atomic
#'           vector or factor of length `nrow(fcs_join_obj$data)`.
#'     \item Named list slots containing `$clusters` or `$coordinates`
#'           (detected automatically).
#'   }
#'   Slots `"pca"` and `"batch_correction"` are intentionally excluded and
#'   will produce a warning if requested (these should be re-computed on the
#'   subsetted object). Non-existent or non-subsettable slots are skipped with
#'   a warning.
#'
#' @details
#'   - If the input object lacks an `object_history` entry, prints a message
#'     recommending `FCSimple::fcs_audit()`.  
#'   - When subsetting by `"source"`, all rows where
#'     `fcs_join_obj$source %in% subset_values` are retained.  
#'   - When subsetting by `"cluster"`, the function extracts
#'     `fcs_join_obj[[subset_cluster_algorithm]]$clusters` and retains rows
#'     matching `subset_values`.  
#'   - The returned object always contains the subsetted `data`, `source`, and
#'     (if present) `run_date` entries. The `raw` slot, if present, is also
#'     automatically subset and retained.  
#'   - A new list element `subset_on` records:
#'     - `subset_by`, `subset_cluster_algorithm`, `subset_values`,
#'       `source_object_indices` (the row indices kept), and
#'       `kept_slots` (names of additional slots that were subset).
#'
#' @return
#'   A list with elements:
#'   - `data`: numeric matrix of selected events × parameters  
#'   - `source`: character vector of sample IDs for selected events  
#'   - `run_date`: character vector of acquisition dates (if original had it)  
#'   - `raw`: numeric matrix of untransformed selected events (if original had it)
#'   - `subset_on`: list recording subsetting criteria, indices, and kept slots  
#'   - Plus any additional slots requested via `keep_slots` (subsetted).
#'
#' @examples
#' \dontrun{
#' # Subset by sample names
#' joined <- FCSimple::fcs_join(files)
#' sel1 <- FCSimple::fcs_subset(
#'   joined,
#'   subset_by = "source",
#'   subset_values = c("SampleA","SampleC")
#' )
#'
#' # Subset by Leiden clusters, keep clustering and UMAP
#' clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#' reduced <- FCSimple::fcs_reduce_dimensions(clustered, algorithm = "umap")
#' sel2 <- FCSimple::fcs_subset(
#'   reduced,
#'   subset_by = "cluster",
#'   subset_cluster_algorithm = "leiden",
#'   subset_values = c(1, 3, 5),
#'   keep_slots = c("leiden", "umap")
#' )
#'
#' # Keep a custom per-cell vector
#' reduced$cell_score <- runif(nrow(reduced$data))
#' sel3 <- FCSimple::fcs_subset(
#'   reduced,
#'   subset_by = "cluster",
#'   subset_cluster_algorithm = "leiden",
#'   subset_values = c(2),
#'   keep_slots = c("leiden", "umap", "cell_score")
#' )
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, FCSimple::fcs_cluster, FCSimple::fcs_audit
#'
#' @export
fcs_subset <- function(fcs_join_obj,
                       subset_by = c("cluster","source"),
                       subset_cluster_algorithm = NA,
                       subset_values,
                       keep_slots = NULL)
{
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_audit() on the object.")
  }

  known_clusterings <- c("leiden", "louvain", "flowsom", "phenograph")
  known_reductions  <- c("umap", "tsne")
  excluded_slots    <- c("pca", "batch_correction")
  n_cells           <- nrow(fcs_join_obj[["data"]])

  if(subset_by=="source") {
    get_index <- which(fcs_join_obj[["source"]] %in% subset_values)
    if(length(get_index)==0) {
      stop(paste0("Value requested not found in source. Use one or more of: ",
                  paste0(names(table(fcs_join_obj[["source"]])), collapse = ", ")))
    }
  } else if(subset_by=="cluster") {
    cluster_nums <- fcs_join_obj[[tolower(subset_cluster_algorithm)]][["clusters"]]
    get_index <- which(cluster_nums %in% subset_values)
  } else {
    stop("error in argument 'subset_by': only subettable by cluster or source")
  }

  # Build the core output
  fcs_new_obj <- list(data   = fcs_join_obj[["data"]][get_index, , drop = FALSE],
                      source = fcs_join_obj[["source"]][get_index])
  if("run_date" %in% names(fcs_join_obj)) {
    fcs_new_obj[["run_date"]] <- fcs_join_obj[["run_date"]][get_index]
  }
  if("raw" %in% names(fcs_join_obj)) {
    if(is.matrix(fcs_join_obj[["raw"]]) || is.data.frame(fcs_join_obj[["raw"]])) {
      fcs_new_obj[["raw"]] <- fcs_join_obj[["raw"]][get_index, , drop = FALSE]
    }
  }

  kept_slots_log <- character(0)

  if(!is.null(keep_slots)) {
    if(!is.character(keep_slots)) {
      stop("'keep_slots' must be a character vector of slot names or NULL")
    }

    for(slot_name in keep_slots) {

      # --- Excluded slots ---
      if(slot_name %in% excluded_slots) {
        warning("Slot '", slot_name, "' is excluded from subsetting. ",
                "Re-run fcs_pca() or fcs_batch_correction() on the subsetted object instead.")
        next
      }

      # --- Slot must exist ---
      if(!slot_name %in% names(fcs_join_obj)) {
        warning("Slot '", slot_name, "' not found in fcs_join_obj. Skipping.")
        next
      }

      slot_content <- fcs_join_obj[[slot_name]]

      # --- 1) Known clustering algorithms ---
      if(slot_name %in% known_clusterings && is.list(slot_content) && "clusters" %in% names(slot_content)) {
        cl_vec <- slot_content[["clusters"]]
        if(length(cl_vec) != n_cells) {
          warning("Slot '", slot_name, "$clusters' length (", length(cl_vec),
                  ") does not match nrow(data) (", n_cells, "). Skipping.")
          next
        }
        fcs_new_obj[[slot_name]] <- slot_content
        fcs_new_obj[[slot_name]][["clusters"]] <- cl_vec[get_index]
        kept_slots_log <- c(kept_slots_log, slot_name)
        next
      }

      # --- 2) Known dimension reductions ---
      if(slot_name %in% known_reductions && is.list(slot_content) && "coordinates" %in% names(slot_content)) {
        coords <- slot_content[["coordinates"]]
        if(!is.matrix(coords) && !is.data.frame(coords)) {
          warning("Slot '", slot_name, "$coordinates' is not a matrix. Skipping.")
          next
        }
        if(nrow(coords) != n_cells) {
          warning("Slot '", slot_name, "$coordinates' has ", nrow(coords),
                  " rows but data has ", n_cells, " rows. Skipping.")
          next
        }
        fcs_new_obj[[slot_name]] <- slot_content
        fcs_new_obj[[slot_name]][["coordinates"]] <- coords[get_index, , drop = FALSE]
        kept_slots_log <- c(kept_slots_log, slot_name)
        next
      }

      # --- 3) Simple per-cell vector ---
      if(is.atomic(slot_content) && !is.matrix(slot_content)) {
        if(length(slot_content) != n_cells) {
          warning("Slot '", slot_name, "' has length ", length(slot_content),
                  " but data has ", n_cells, " rows. Skipping.")
          next
        }
        fcs_new_obj[[slot_name]] <- slot_content[get_index]
        kept_slots_log <- c(kept_slots_log, slot_name)
        next
      }

      # --- 4) Named list slot with $clusters ---
      if(is.list(slot_content) && "clusters" %in% names(slot_content)) {
        cl_vec <- slot_content[["clusters"]]
        if(length(cl_vec) != n_cells) {
          warning("Slot '", slot_name, "$clusters' length (", length(cl_vec),
                  ") does not match nrow(data) (", n_cells, "). Skipping.")
          next
        }
        fcs_new_obj[[slot_name]] <- slot_content
        fcs_new_obj[[slot_name]][["clusters"]] <- cl_vec[get_index]
        kept_slots_log <- c(kept_slots_log, slot_name)
        next
      }

      # --- 5) Named list slot with $coordinates (matrix) ---
      if(is.list(slot_content) && "coordinates" %in% names(slot_content)) {
        coords <- slot_content[["coordinates"]]
        if(is.matrix(coords) || is.data.frame(coords)) {
          if(nrow(coords) != n_cells) {
            warning("Slot '", slot_name, "$coordinates' has ", nrow(coords),
                    " rows but data has ", n_cells, " rows. Skipping.")
            next
          }
          fcs_new_obj[[slot_name]] <- slot_content
          fcs_new_obj[[slot_name]][["coordinates"]] <- coords[get_index, , drop = FALSE]
          kept_slots_log <- c(kept_slots_log, slot_name)
          next
        }
      }

      # --- 6) Cannot subset ---
      warning("Slot '", slot_name, "' could not be subset (unrecognized structure). Skipping.")
    }
  }

  fcs_new_obj[["subset_on"]] <- list(subset_by = subset_by,
                                     subset_cluster_algorithm = subset_cluster_algorithm,
                                     subset_values = subset_values,
                                     source_object_indices = get_index,
                                     kept_slots = kept_slots_log)
  return(fcs_new_obj)
}

#' @title Remove Specified Parameters from FCSimple Object
#'
#' @description
#'   Drops one or more channels (columns) from the `data` matrix of an
#'   FCSimple object, returning a new object with the reduced parameter set.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), containing at least a
#'   `data` matrix and a `source` vector. May also include `run_date`.
#'
#' @param remove_parameters
#'   Character vector of column names (parameters) to remove from
#'   `fcs_join_obj$data`.
#'
#' @details
#'   - Searches `colnames(fcs_join_obj$data)` for any entries in
#'     `remove_parameters`.  
#'   - If none are found, the function errors.  
#'   - Otherwise, all matching columns are dropped.
#'
#' @return
#'   A list with elements:
#'   - `data`: numeric matrix of events × remaining parameters  
#'   - `source`: character vector of sample IDs (unchanged)  
#'   - `run_date`: character vector of acquisition dates (if original had it)  
#'
#' @examples
#' \dontrun{
#' joined <- FCSimple::fcs_join(files)
#' # Remove CD3 and CD19 channels before downstream analysis
#' pruned <- FCSimple::fcs_remove_parameters(
#'   joined,
#'   remove_parameters = c("CD3","CD19")
#' )
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, FCSimple::fcs_subset
#'
#' @export
fcs_remove_parameters <- function(fcs_join_obj,
                                  remove_parameters)
{
  rm_param <- which(colnames(fcs_join_obj[["data"]]) %in% remove_parameters)
    if(length(rm_param) > 0) {
      fcs_join_obj$data <- fcs_join_obj$data[,-rm_param]
      if('raw' %in% names(fcs_join_obj)) {
        fcs_join_obj$raw <- fcs_join_obj$raw[,-rm_param]
      }
    } else {
      message("No parameters were removed.")
    }
    return(fcs_join_obj)
}
