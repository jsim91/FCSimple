#' @title Quantify Semi-Local Neighborhood Consistency Across Reductions
#'
#' @description
#'   Compares each low-dimensional embedding (UMAP/t-SNE, 2D/3D) against the
#'   high-dimensional k-nearest-neighbor graph used for clustering.  For every
#'   cell it measures how faithfully the low-dimensional embedding preserves
#'   the true high-dimensional semi-local neighborhood — where "neighborhood"
#'   is an \emph{r-hop} set grown by graph traversal over the k-NN graphs
#'   computed by \code{FCSimple::fcs_embedding_neighbors()}.
#'
#'   The r-hop form makes the metric robust to dataset size: on millions of
#'   cells a fixed k captures only an infinitesimal patch of the manifold, so
#'   any projection wobble scrambles that microscopic set and all scores
#'   collapse toward zero.  Expanding to r hops (e.g. 2-3) yields semi-local
#'   neighborhoods that stay stable under subtle landscape shifts while still
#'   detecting gross topological breaks.
#'
#' @param fcs_join_obj
#'   A list whose \code{$neighbors} element (from
#'   \code{FCSimple::fcs_embedding_neighbors()}) includes \code{"high_dim"}
#'   and one entry per embedding, each with a \code{hop_idx} list of r-hop
#'   neighbor cell indices.
#'
#' @param n_hops
#'   Integer >= 1; the hop radius to use (must have been computed during
#'   \code{fcs_embedding_neighbors()}).  Default \code{1}.
#'
#' @details
#'   The high-dimensional reference is \code{fcs_join_obj$neighbors$high_dim},
#'   and each embedding's set is \code{fcs_join_obj$neighbors[[slot]]$hop_idx}.
#'   Because both are produced by the identical graph-traversal operation with
#'   the same \code{n_hops}, the sets are definition- and cardinality-matched.
#'
#'   Three complementary scores are computed per cell:
#'   \describe{
#'     \item{\code{jaccard}}{Intersection size divided by union size of the
#'       high- and low-dimensional r-hop sets (0-1).}
#'     \item{\code{continuity}}{Intersection size divided by the high-dimensional
#'       r-hop set size: the fraction of true high-dimensional neighbours kept
#'       close.  Low = neighbours \emph{torn apart}.}
#'     \item{\code{trustworthiness}}{Intersection size divided by the
#'       low-dimensional r-hop set size: the fraction of low-dimensional
#'       neighbours that are genuine.  Low = \emph{false adjacencies}.}
#'   }
#'   All scores lie in \eqn{[0, 1]} with 1 denoting perfect preservation.
#'
#' @return
#'   The input \code{fcs_join_obj} augmented with
#'   \code{$neighborhood_consistency}, a named list with one entry per
#'   embedding (e.g. \code{"umap_2d"}, \code{"umap_3d"}, \code{"tsne_2d"},
#'   \code{"tsne_3d"}).  Each entry is a data.frame with columns
#'   \code{jaccard}, \code{continuity}, \code{trustworthiness}, and
#'   \code{nn_mean_dist}.
#'
#'   \code{object_history} is appended with a timestamped entry.
#'
#' @examples
#' \dontrun{
#'   obj <- fcs_embedding_neighbors(obj, k = 30, n_hops = 2)
#'   obj <- fcs_neighborhood_consistency(obj, n_hops = 2)
#'   hist(obj$neighborhood_consistency$umap_2d$trustworthiness)
#' }
#'
#' @seealso
#'   FCSimple::fcs_embedding_neighbors, FCSimple::fcs_plot_consistency
#'
#' @export
fcs_neighborhood_consistency <- function(fcs_join_obj, n_hops = 1L) {
  neighbors <- fcs_join_obj$neighbors
  if (is.null(neighbors))
    stop("No '$neighbors' found. Run fcs_embedding_neighbors() first.")
  if (!"high_dim" %in% names(neighbors) || is.null(neighbors$high_dim))
    stop("No 'high_dim' entry in $neighbors. Run fcs_embedding_neighbors() ",
         "after fcs_cluster().")

  n_hops <- as.integer(n_hops)
  if (n_hops < 1L) stop("'n_hops' must be >= 1.")

  n_cells <- length(neighbors$high_dim$nn_mean_dist)
  high_hop <- neighbors$high_dim$hop_idx
  if (is.null(high_hop)) {
    stop("No 'hop_idx' in $neighbors$high_dim. Re-run ",
         "fcs_embedding_neighbors() with the desired n_hops.")
  }

  emb_slots <- setdiff(names(neighbors), "high_dim")
  out <- list()

  for (slot in emb_slots) {
    low_hop <- neighbors[[slot]]$hop_idx
    if (is.null(low_hop)) next

    inter  <- vapply(seq_len(n_cells), function(i) {
      length(intersect(high_hop[[i]], low_hop[[i]]))
    }, integer(1))
    hi_sz  <- lengths(high_hop)
    lo_sz  <- lengths(low_hop)

    union_sz        <- hi_sz + lo_sz - inter
    jaccard         <- ifelse(union_sz > 0L, inter / union_sz, 0)
    continuity      <- ifelse(hi_sz > 0L,    inter / hi_sz,    0)
    trustworthiness <- inter / pmax(lo_sz, 1L)

    out[[slot]] <- data.frame(
      jaccard         = jaccard,
      continuity      = continuity,
      trustworthiness = trustworthiness,
      nn_mean_dist    = neighbors[[slot]]$nn_mean_dist,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(fcs_join_obj$neighborhood_consistency)) {
    fcs_join_obj$neighborhood_consistency <- out
  } else {
    for (s in names(out)) {
      fcs_join_obj$neighborhood_consistency[[s]] <- out[[s]]
    }
  }

  if (!"object_history" %in% names(fcs_join_obj)) {
    message("Consider running FCSimple::fcs_audit() on the object.")
  }
  try(expr = fcs_join_obj[["object_history"]] <- append(
    fcs_join_obj[["object_history"]],
    paste0("neighborhood_consistency (n_hops=", n_hops, "; ",
           paste(names(out), collapse = ", "), "): ", Sys.time())
  ), silent = TRUE)

  message("Neighbourhood consistency computed for ", length(out),
          " embedding(s) at n_hops = ", n_hops, ".")
  return(fcs_join_obj)
}