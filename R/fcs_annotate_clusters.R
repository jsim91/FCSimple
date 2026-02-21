#' @title Annotate Clusters with Cell Type Labels
#'
#' @description
#'   Assigns cell type names to clusters in an FCSimple analysis object by
#'   building a `{algorithm}_mapping` data frame (e.g. `leiden_mapping`,
#'   `flowsom_mapping`). This mapping is used by the FCView Shiny application
#'   to populate the Annotation tab and enable cell-type-level analysis across
#'   all downstream tabs.
#'
#' @param fcs_join_obj
#'   An FCSimple analysis object (list) containing at minimum a clustering
#'   result with a `clusters` vector (e.g. from `FCSimple::fcs_cluster()`).
#'   The clustering result must be stored under a recognized algorithm name
#'   (`"leiden"`, `"flowsom"`, `"louvain"`, `"phenograph"`) or `"cluster"`.
#'
#' @param annotations
#'   A named list mapping cell type labels to one or more cluster IDs.
#'   Names are the cell type labels; values are integer or character vectors
#'   of cluster IDs to assign to that label. Example:
#'   \preformatted{
#'   list(
#'     "T cell" = c(1, 2, 3),
#'     "B cell" = c(4, 6, 7),
#'     "NK"     = 5
#'   )
#'   }
#'   Not all clusters need to be annotated. Any cluster not referenced will
#'   receive a `NA` value in the `celltype` column of the mapping, indicating
#'   it is unassigned. The FCView annotation tab treats `NA`-mapped clusters
#'   as pending assignment and will not pre-fill a cell type for them.
#'
#' @param clustering_algorithm
#'   Character or `NULL` (default). Name of the clustering algorithm whose
#'   cluster IDs are referenced in `annotations`
#'   (e.g. `"leiden"`, `"flowsom"`, `"cluster"`). If `NULL`, the first
#'   detected algorithm is used. The resulting mapping is stored on the
#'   object under this name: e.g. `"leiden"` → `leiden_mapping`. This
#'   allows multiple algorithm mappings to coexist on the same object.
#'
#' @param overwrite
#'   Logical; if `TRUE` (default), replaces any existing `{algorithm}_mapping`
#'   element for the selected algorithm. If `FALSE` and the mapping already
#'   exists, an error is raised.
#'
#' @details
#'   The function:
#'   \enumerate{
#'     \item Validates that `annotations` is a non-empty named list.
#'     \item Detects or validates the specified clustering algorithm.
#'     \item Determines the set of valid cluster IDs from the cluster vector.
#'     \item Warns about any annotated cluster IDs not present in the data.
#'     \item Warns about any cluster IDs assigned to multiple cell types
#'           (each cluster can only belong to one cell type; the first
#'           assignment wins).
#'     \item Builds a two-column data frame with columns `cluster` (integer)
#'           and `celltype` (character), one row per cluster in the data,
#'           sorted by cluster ID.
#'     \item Clusters not present in `annotations` receive `NA` in the
#'           `celltype` column, flagging them as unassigned.
#'     \item Stores the result as `fcs_join_obj${algorithm}_mapping`
#'           (e.g. `leiden_mapping`, `flowsom_mapping`).
#'   }
#'
#'   The FCView annotation tab will not pre-fill a cell type for `NA`-mapped
#'   clusters, allowing the user to assign them interactively in the app.
#'
#' @return
#'   The input `fcs_join_obj` with a `cluster_mapping` data frame added (or
#'   replaced), and an updated `object_history` entry.
#'
#' @examples
#' \dontrun{
#'   files   <- FCSimple::fcs_example_files()
#'   joined  <- FCSimple::fcs_join(files)
#'   clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#'
#'   # Annotate a subset of clusters (others left un-named)
#'   annotated <- FCSimple::fcs_annotate_clusters(
#'     clustered,
#'     annotations = list(
#'       "CD4 T cell" = c(1, 3),
#'       "CD8 T cell" = c(2, 5),
#'       "B cell"     = 4,
#'       "NK"         = 6
#'     )
#'   )
#'
#'   annotated$cluster_mapping
#'   #   cluster    celltype
#'   # 1       1  CD4 T cell
#'   # 2       2  CD8 T cell
#'   # 3       3  CD4 T cell
#'   # 4       4      B cell
#'   # 5       5  CD8 T cell
#'   # 6       6          NK
#'
#'   # Prepare for FCView (cluster_mapping is preserved automatically)
#'   # The mapping is stored as leiden_mapping because clustering_algorithm
#'   # defaults to the first detected algorithm (leiden here).
#'   names(annotated)  # includes "leiden_mapping"
#'
#'   # Multiple algorithms: annotate each separately
#'   clustered2 <- FCSimple::fcs_cluster(clustered, algorithm = "flowsom")
#'   annotated2 <- FCSimple::fcs_annotate_clusters(
#'     clustered2,
#'     annotations = list("T cell" = c(1, 2), "B cell" = 3),
#'     clustering_algorithm = "flowsom"
#'   )
#'   # annotated2$flowsom_mapping is now available alongside leiden_mapping
#'
#'   # Prepare for FCView — pass the algorithm whose mapping you want to use:
#'   prepared <- FCSimple::fcs_prepare_fcview_object(
#'     annotated,
#'     clustering_algorithm = "leiden"
#'   )
#'   # prepared$cluster_mapping is leiden_mapping renamed for the app
#' }
#'
#' @seealso
#'   FCSimple::fcs_cluster, FCSimple::fcs_prepare_fcview_object
#'
#' @export
fcs_annotate_clusters <- function(fcs_join_obj,
                                  annotations,
                                  clustering_algorithm = NULL,
                                  overwrite = TRUE) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!is.list(fcs_join_obj)) {
    stop("fcs_join_obj must be a list")
  }

  if (!is.list(annotations) || length(annotations) == 0) {
    stop("annotations must be a non-empty named list")
  }

  ann_names <- names(annotations)
  if (is.null(ann_names) || any(!nzchar(ann_names))) {
    stop("All elements of annotations must be named (non-empty character names)")
  }

  # overwrite check is deferred until selected_algo is known (see below)

  # ── Detect clustering algorithm ─────────────────────────────────────────────
  all_cluster_algos <- c("leiden", "flowsom", "louvain", "phenograph", "cluster")
  present_algos <- Filter(function(algo) {
    algo %in% names(fcs_join_obj) &&
      is.list(fcs_join_obj[[algo]]) &&
      "clusters" %in% names(fcs_join_obj[[algo]])
  }, all_cluster_algos)

  if (length(present_algos) == 0) {
    stop("No clustering result found in fcs_join_obj. Expected one of: ",
         paste(all_cluster_algos, collapse = ", "))
  }

  if (!is.null(clustering_algorithm)) {
    if (!clustering_algorithm %in% present_algos) {
      stop("Specified clustering_algorithm '", clustering_algorithm,
           "' not found. Available: ", paste(present_algos, collapse = ", "))
    }
    selected_algo <- clustering_algorithm
  } else {
    selected_algo <- present_algos[1]
    if (length(present_algos) > 1) {
      message("Multiple clustering algorithms detected: ",
              paste(present_algos, collapse = ", "))
      message("Using '", selected_algo, "'. ",
              "Specify clustering_algorithm to choose a different one.")
    }
  }

  mapping_name   <- paste0(selected_algo, "_mapping")
  if (!overwrite && !is.null(fcs_join_obj[[mapping_name]])) {
    stop(mapping_name, " already exists. Set overwrite = TRUE to replace it.")
  }

  cluster_vec    <- fcs_join_obj[[selected_algo]]$clusters
  valid_clusters <- unique(cluster_vec)

  # ── Build cluster → celltype mapping ────────────────────────────────────────
  # Walk annotations once; first assignment of a cluster wins (with a warning
  # if a cluster appears in more than one cell type).
  seen_clusters <- integer(0)
  rows          <- list()

  for (celltype in ann_names) {
    ids <- as.integer(annotations[[celltype]])

    # Warn about IDs not present in the cluster vector
    bad_ids <- setdiff(ids, valid_clusters)
    if (length(bad_ids) > 0) {
      warning(sprintf(
        "Cell type '%s': cluster ID(s) not found in data and will be skipped: %s",
        celltype, paste(bad_ids, collapse = ", ")
      ))
      ids <- intersect(ids, valid_clusters)
    }

    if (length(ids) == 0) next

    # Warn about clusters already assigned to a different cell type
    dupes <- intersect(ids, seen_clusters)
    if (length(dupes) > 0) {
      warning(sprintf(
        "Cluster ID(s) %s already assigned to an earlier cell type; skipping duplicate assignment to '%s'.",
        paste(dupes, collapse = ", "), celltype
      ))
      ids <- setdiff(ids, dupes)
    }

    if (length(ids) == 0) next

    seen_clusters <- c(seen_clusters, ids)
    rows <- c(rows, list(data.frame(
      cluster  = ids,
      celltype = celltype,
      stringsAsFactors = FALSE
    )))
  }

  if (length(rows) == 0) {
    stop("No valid cluster annotations after validation. Check that cluster IDs match those in the object.")
  }

  cluster_mapping <- do.call(rbind, rows)
  cluster_mapping <- cluster_mapping[order(cluster_mapping$cluster), ]
  rownames(cluster_mapping) <- NULL

  # Add NA rows for every cluster not covered by annotations
  all_cluster_ids  <- as.integer(sort(unique(valid_clusters)))
  unannotated_ids  <- setdiff(all_cluster_ids, cluster_mapping$cluster)
  n_unannotated    <- length(unannotated_ids)

  if (n_unannotated > 0) {
    na_rows <- data.frame(
      cluster  = unannotated_ids,
      celltype = NA_character_,
      stringsAsFactors = FALSE
    )
    cluster_mapping <- rbind(cluster_mapping, na_rows)
    cluster_mapping <- cluster_mapping[order(cluster_mapping$cluster), ]
    rownames(cluster_mapping) <- NULL
  }

  # ── Store result as {algo}_mapping ───────────────────────────────────────────
  fcs_join_obj[[mapping_name]] <- cluster_mapping

  n_annotated <- sum(!is.na(cluster_mapping$celltype))
  n_total     <- nrow(cluster_mapping)

  message(sprintf(
    "%s created: %d/%d cluster(s) annotated as %d cell type(s) using algorithm '%s'.",
    mapping_name, n_annotated, n_total,
    length(unique(stats::na.omit(cluster_mapping$celltype))), selected_algo
  ))
  if (n_unannotated > 0) {
    message(sprintf(
      "  %d cluster(s) left un-annotated (celltype = NA; can be assigned in FCView).",
      n_unannotated
    ))
  }

  # ── History ──────────────────────────────────────────────────────────────────
  fcs_join_obj$object_history <- c(
    fcs_join_obj$object_history,
    list(list(
      step                 = "fcs_annotate_clusters",
      timestamp            = Sys.time(),
      algorithm            = selected_algo,
      mapping_element      = mapping_name,
      n_celltypes          = length(unique(stats::na.omit(cluster_mapping$celltype))),
      n_annotated_clusters = n_annotated,
      n_total_clusters     = n_total
    ))
  )

  return(fcs_join_obj)
}
