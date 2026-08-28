#' @title Meta-Region Composition Analysis Across Reductions
#'
#' @description
#'   Partitions the high-dimensional k-NN graph into many small "meta-regions"
#'   using high-resolution Leiden clustering, then measures how faithfully each
#'   meta-region's composition is preserved when cells are projected into a
#'   low-dimensional embedding.  This provides a coarser, more interpretable
#'   view of high-to-low-dimensional structural integrity than per-cell
#'   neighbourhood consistency alone.
#'
#' @param fcs_join_obj
#'   A list produced by \code{FCSimple::fcs_reduce_dimensions()},
#'   \code{FCSimple::fcs_cluster()}, and
#'   \code{FCSimple::fcs_embedding_neighbors()}.  The high-dimensional
#'   adjacency (or search) matrix must be present.
#'
#' @param resolution_parameter
#'   Numeric; Leiden resolution for the meta-region partition.  Higher values
#'   produce more, smaller meta-regions.  Default \code{5}.
#'
#' @param n_hops
#'   Integer >= 1; hop radius of the semi-local neighbourhood used to assess
#'   low-dimensional mixing.  Must match a radius already computed by
#'   \code{fcs_embedding_neighbors()}.  Default \code{1}.
#'
#' @param seed
#'   Integer; random seed for Leiden (where applicable).  Default \code{123}.
#'
#' @details
#'   \enumerate{
#'     \item The high-dimensional graph is rebuilt from
#'           \code{fcs_join_obj$adjacency_matrix} (or \code{$search}) and
#'           partitioned with high-resolution \code{igraph::cluster_leiden()}.
#'     \item For each resulting meta-region \code{m}, and for each low-dimensional
#'           embedding, two integrity scores are computed over the member cells'
#'           r-hop neighbourhoods:
#'       \describe{
#'         \item{\code{retention}}{Mean fraction of a member cell's low-D r-hop
#'           neighbours that belong to the \emph{same} meta-region.  High =
#'           the region stays coherent in low dimensions.}
#'         \item{\code{mixing}}{Mean fraction of foreign cells (from other
#'           meta-regions) invading a member's low-D r-hop neighbourhood.
#'           High = the region is dispersed in low dimensions.}
#'       }
#'     \item A per-embedding summary of mean retention/mixing is also stored.
#'   }
#'
#'   These scores are \emph{informative at the semi-local scale}, complementing
#'   \code{fcs_neighborhood_consistency()}: while the latter flags individual
#'   torn or invented edges, meta-region retention/mixing reports whether whole
#'   coherent high-dimensional patches survive the projection.
#'
#' @return
#'   The input \code{fcs_join_obj} augmented with \code{$meta_regions}, a list:
#'   \describe{
#'     \item{\code{membership}}{Integer vector; meta-region id per cell.}
#'     \item{\code{regions}}{Named list (by embedding slot) of per-meta-region
#'       data.frames with columns \code{meta_region}, \code{n_cells},
#'       \code{retention}, and \code{mixing}.}
#'     \item{\code{summary}}{data.frame with one row per embedding summarising
#'       mean retention and mixing across meta-regions.}
#'     \item{\code{features}}{Character vector of the features the partitioned
#'       high-dimensional graph was built from (\code{$search_features} or
#'       \code{$neighbors$high_dim$features}).}
#'   }
#'
#'   \code{object_history} is appended with a timestamped entry.
#'
#' @examples
#' \dontrun{
#'   obj <- fcs_embedding_neighbors(obj, k = 30, n_hops = 2)
#'   obj <- fcs_meta_regions(obj, resolution_parameter = 5, n_hops = 2)
#'   head(obj$meta_regions$summary)
#' }
#'
#' @seealso
#'   FCSimple::fcs_embedding_neighbors, FCSimple::fcs_neighborhood_consistency
#'
#' @importFrom igraph graph.adjacency cluster_leiden make_empty_graph add_edges
#' @importFrom Matrix summary
#' @export
fcs_meta_regions <- function(
    fcs_join_obj,
    resolution_parameter = 5,
    n_hops = 1L,
    seed = 123L)
{
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required. Install with: install.packages('igraph')")

  neighbors <- fcs_join_obj$neighbors
  if (is.null(neighbors))
    stop("No '$neighbors' found. Run fcs_embedding_neighbors() first.")
  if (!"high_dim" %in% names(neighbors) || is.null(neighbors$high_dim))
    stop("No 'high_dim' entry in $neighbors. Run fcs_embedding_neighbors() ",
         "after fcs_cluster().")

  n_hops <- as.integer(n_hops)
  if (n_hops < 1L) stop("'n_hops' must be >= 1.")

  hd <- neighbors$high_dim
  n_cells <- length(hd$nn_mean_dist)

  # -- Build high-dimensional undirected graph --------------------------------
  if (!is.null(fcs_join_obj$adjacency_matrix)) {
    coords <- Matrix::summary(fcs_join_obj$adjacency_matrix)
    ii <- coords$i
    jj <- coords$j
  } else if (!is.null(hd$nn_i)) {
    ii <- hd$nn_i
    jj <- hd$nn_j
  } else if (!is.null(hd$nn_idx)) {
    kh <- ncol(hd$nn_idx)
    ii <- rep(seq_len(n_cells), each = kh)
    jj <- as.vector(t(hd$nn_idx))
  } else {
    stop("No high-dimensional graph representation available.")
  }

  # Build the undirected graph from the edge list (robust and memory-light).
  # The adjacency/search matrix is already an undirected, non-negative graph
  # from fcs_cluster(), so the edge list is used directly.
  g <- igraph::make_empty_graph(n = n_cells, directed = FALSE)
  g <- igraph::add_edges(g, edges = as.vector(rbind(ii, jj)))

  set.seed(seed)
  message("Running high-resolution Leiden (resolution = ",
          resolution_parameter, ")...")
  part <- igraph::cluster_leiden(
    g,
    objective_function = "modularity",
    resolution_parameter = resolution_parameter,
    weights = NA
  )
  membership <- as.integer(part$membership)
  meta_ids <- sort(unique(membership))
  n_meta <- length(meta_ids)

  message("Meta-region partition complete: ", n_meta, " meta-regions for ",
          n_cells, " cells.")

  # -- Low-dimensional integrity per embedding --------------------------------
  emb_slots <- setdiff(names(neighbors), "high_dim")
  regions <- list()
  summary <- list()

  for (slot in emb_slots) {
    low_hop <- neighbors[[slot]]$hop_idx
    if (is.null(low_hop)) next

    region_stats <- lapply(meta_ids, function(mid) {
      cells <- which(membership == mid)
      if (length(cells) == 0L)
        return(data.frame(meta_region = mid, n_cells = 0L,
                          retention = NA_real_, mixing = NA_real_))
      # For each member cell, compute same-meta vs foreign fraction in its
      # low-D r-hop neighborhood.
      ret <- vapply(cells, function(ci) {
        h <- low_hop[[ci]]
        if (length(h) == 0L) return(NA_real_)
        mean(membership[h] == mid)
      }, numeric(1))
      mix <- vapply(cells, function(ci) {
        h <- low_hop[[ci]]
        if (length(h) == 0L) return(NA_real_)
        mean(membership[h] != mid)
      }, numeric(1))
      data.frame(meta_region = mid, n_cells = length(cells),
                 retention = mean(ret, na.rm = TRUE),
                 mixing    = mean(mix, na.rm = TRUE))
    })
    regions[[slot]] <- do.call(rbind, region_stats)

    summary[[slot]] <- data.frame(
      slot            = slot,
      n_meta_regions  = n_meta,
      mean_retention  = mean(regions[[slot]]$retention,  na.rm = TRUE),
      mean_mixing     = mean(regions[[slot]]$mixing,     na.rm = TRUE)
    )
  }

  summary_df <- do.call(rbind, summary)
  rownames(summary_df) <- NULL

  meta_features <- fcs_join_obj$search_features
  if (is.null(meta_features) && !is.null(neighbors$high_dim$features)) {
    meta_features <- neighbors$high_dim$features
  }
  fcs_join_obj$meta_regions <- list(
    membership = membership,
    regions    = regions,
    summary    = summary_df,
    features   = meta_features,
    resolution_parameter = resolution_parameter,
    n_hops     = n_hops
  )

  if (!"object_history" %in% names(fcs_join_obj)) {
    message("Consider running FCSimple::fcs_audit() on the object.")
  }
  try(expr = fcs_join_obj[["object_history"]] <- append(
    fcs_join_obj[["object_history"]],
    paste0("meta_regions (res=", resolution_parameter, ", n_hops=", n_hops,
           "): ", Sys.time())
  ), silent = TRUE)

  cat("\nMeta-region integrity summary:\n")
  print(summary_df)
  return(fcs_join_obj)
}