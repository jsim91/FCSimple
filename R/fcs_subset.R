#' @title Subset FCSimple Analysis Object by Sample or Cluster
#'
#' @description
#'   Creates a new FCSimple-style object containing only the events that match
#'   specified sample names or cluster IDs. Records the subsetting criteria and
#'   original indices in a `subset_on` element for provenance.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), optionally processed by
#'   FCSimple::fcs_cluster() or other functions. Should include at least
#'   `data` (numeric matrix) and `source` (vector of sample identifiers).
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
#' @details
#'   - If the input object lacks an `object_history` entry, prints a message
#'     recommending `FCSimple::fcs_audit()`.  
#'   - When subsetting by `"source"`, all rows where
#'     `fcs_join_obj$source %in% subset_values` are retained.  
#'   - When subsetting by `"cluster"`, the function extracts
#'     `fcs_join_obj[[subset_cluster_algorithm]]$clusters` and retains rows
#'     matching `subset_values`.  
#'   - The returned object contains only the selected `data`, `source`, and
#'     (if present) `run_date` entries.  
#'   - A new list element `subset_on` records:
#'     - `subset_by`, `subset_cluster_algorithm`, `subset_values`, and
#'       `source_object_indices` (the row indices kept).
#'
#' @return
#'   A list with elements:
#'   - `data`: numeric matrix of selected events × parameters  
#'   - `source`: character vector of sample IDs for selected events  
#'   - `run_date`: character vector of acquisition dates (if original had it)  
#'   - `subset_on`: list recording subsetting criteria and indices  
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
#' # Subset by Leiden clusters
#' clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#' sel2 <- FCSimple::fcs_subset(
#'   clustered,
#'   subset_by = "cluster",
#'   subset_cluster_algorithm = "leiden",
#'   subset_values = c(1, 3, 5)
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
                       subset_values)
{
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_audit() on the object.")
  }
  if(subset_by=="source") {
    get_index <- which(fcs_join_obj[["source"]] %in% subset_values)
    if(length(get_index)==0) {
      stop(paste0("Value requested not found in source. Use one or more of: ",
                  paste0(names(table(fcs_join_obj[["source"]])), collapse = ", ")))
    }
    fcs_new_obj <- list(data = fcs_join_obj[["data"]][get_index,],
                        source = fcs_join_obj[["source"]][get_index])
    if("run_date" %in% names(fcs_join_obj)) {
      fcs_new_obj[["run_date"]] <- fcs_join_obj[["run_date"]][get_index]
    }
  } else if(subset_by=="cluster") {
    cluster_nums <- fcs_join_obj[[tolower(subset_cluster_algorithm)]][["clusters"]]
    get_index <- which(cluster_nums %in% subset_values)
    fcs_new_obj <- list(data = fcs_join_obj[["data"]][get_index,],
                        source = fcs_join_obj[["source"]][get_index])
    if("run_date" %in% names(fcs_join_obj)) {
      fcs_new_obj[["run_date"]] <- fcs_join_obj[["run_date"]][get_index]
    }
  } else {
    stop("error in argument 'subset_by': only subettable by cluster or source")
  }
  fcs_new_obj[["subset_on"]] <- list(subset_by = subset_by,
                                     subset_cluster_algorithm = subset_cluster_algorithm,
                                     subset_values = subset_values,
                                     source_object_indices = get_index)
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
