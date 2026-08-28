#' @title Detect Border Points in Reduction Embeddings
#'
#' @description
#'   Computes border scores for each reduction embedding using the
#'   pre-computed local neighbor graphs stored in \code{fcs_join_obj$neighbors}
#'   (produced by \code{FCSimple::fcs_embedding_neighbors()}).  Two
#'   complementary measures are derived: an inter-cluster border score
#'   (\code{frac_diff}) and an outer-periphery / skin score
#'   (\code{periphery_score}), combined into a single \code{border_score}.
#'
#' @param fcs_join_obj
#'   A list produced by \code{FCSimple::fcs_cluster()} and
#'   \code{FCSimple::fcs_embedding_neighbors()}.  Must contain a
#'   \code{$neighbors} element with one entry per embedding and, optionally, a
#'   \code{"high_dim"} entry.  Cluster labels are taken from the clustering
#'   result named by \code{algorithm}.
#'
#' @param algorithm
#'   Character; name of the clustering result to use for border detection
#'   (e.g. \code{"leiden"}, \code{"flowsom"}).
#'
#' @param n_hops
#'   Integer >= 1; hop radius used for the semi-local \code{frac_diff} border
#'   score.  Must match a radius already computed by
#'   \code{FCSimple::fcs_embedding_neighbors()}.  Default \code{1}.
#'
#' @param num_cores
#'   Integer; kept for API parity with the neighbour computation but unused
#'   here (scoring is single-threaded).  Default
#'   \code{ceiling(parallel::detectCores() / 2)}.
#'
#' @details
#'   Border scoring is purely derived from the stored neighbour indices and
#'   distances; no neighbour search is performed by this function.
#'   \describe{
#'     \item{\code{frac_diff}}{
#'       Fraction of a cell's r-hop semi-local neighbourhood belonging to a
#'       different cluster.  Core cells ≈ 0; two-cluster boundary ≈ 0.5;
#'       three-cluster junction ≈ 0.67.}
#'     \item{\code{periphery_score}}{
#'       Outer-edge score of the full point cloud, independent of cluster
#'       membership.  In 2D this is the largest empty angular sector around
#'       the cell (normalised by \eqn{2\pi}); in 3D it is the resultant norm
#'       of unit vectors to all neighbours.  Interior cells ≈ 0;
#'       outer-envelope cells ≈ 1.}
#'     \item{\code{border_score}}{
#'       \code{pmax(frac_diff, periphery_score)}.}
#'   }
#'   For the high-dimensional reference, only \code{frac_diff} is available
#'   (\code{periphery_score} is \code{NA}).
#'
#' @return
#'   The input \code{fcs_join_obj} augmented with \code{$reduction_borders}, a
#'   named list with one entry per embedding (e.g. \code{"umap_2d"},
#'   \code{"umap_3d"}, \code{"tsne_2d"}, \code{"tsne_3d"}) plus a
#'   \code{"high_dim"} entry when clustering was performed.  Each entry is a
#'   data.frame with columns:
#'   \describe{
#'     \item{\code{frac_diff}}{Numeric; inter-cluster border score (0-1).}
#'     \item{\code{periphery_score}}{Numeric; outer-periphery score (0-1).}
#'     \item{\code{border_score}}{Numeric; combined score.}
#'     \item{\code{nn_mean_dist}}{Numeric; mean neighbour distance.}
#'     \item{\code{cluster}}{Character; cluster ID per cell.}
#'   }
#'
#'   \code{object_history} is appended with a timestamped entry.
#'
#' @examples
#' \dontrun{
#'   obj <- FCSimple::fcs_embedding_neighbors(obj, k = 30)
#'   obj <- FCSimple::fcs_reduction_borders(obj, "leiden")
#'   head(obj$reduction_borders$umap_2d)
#' }
#'
#' @seealso
#'   FCSimple::fcs_embedding_neighbors, FCSimple::fcs_plot_borders
#'
#' @export
fcs_reduction_borders <- function(
    fcs_join_obj,
    algorithm,
    n_hops = 1L,
    num_cores = ceiling(parallel::detectCores() / 2))
{
  # -- Validate clustering ---------------------------------------------------
  if (!tolower(algorithm) %in% names(fcs_join_obj))
    stop("Clustering '", algorithm, "' not found in fcs_join_obj.")

  cluster_labels <- fcs_join_obj[[tolower(algorithm)]][["clusters"]]
  if (is.null(cluster_labels))
    stop("No cluster labels found under fcs_join_obj$", tolower(algorithm),
         "$clusters.")
  cluster_labels <- as.character(cluster_labels)
  n_cells <- length(cluster_labels)

  # -- Require pre-computed neighbors ----------------------------------------
  neighbors <- fcs_join_obj$neighbors
  if (is.null(neighbors))
    stop("No '$neighbors' found. Run fcs_embedding_neighbors() first.")

  emb_slots <- setdiff(names(neighbors), "high_dim")
  if (length(emb_slots) == 0L)
    stop("No embedding neighbour graphs found in $neighbors.")

  border_results <- list()

  # -- Score each embedding --------------------------------------------------
  for (slot in emb_slots) {
    nb <- neighbors[[slot]]
    coords <- fcs_join_obj[[slot]]$coordinates
    if (!is.null(coords) && !is.matrix(coords))
      coords <- as.matrix(coords)

    nn_idx   <- nb$nn_idx
    nn_dist  <- nb$nn_dist
    n_dims   <- if (!is.null(coords)) ncol(coords) else ncol(nn_idx)
    k        <- ncol(nn_idx)

    # frac_diff: fraction of the r-hop semi-local neighbourhood in a
    # different cluster.  Uses the pre-computed hop sets so the border score
    # can be tuned to a semi-local scale independent of k.
    hop_idx <- nb$hop_idx
    if (!is.null(hop_idx) && length(hop_idx) == n_cells) {
      frac_diff <- vapply(seq_len(n_cells), function(i) {
        h <- hop_idx[[i]]
        if (length(h) == 0L) return(0)
        mean(cluster_labels[h] != cluster_labels[i])
      }, numeric(1))
    } else {
      nn_clusters <- matrix(cluster_labels[nn_idx], nrow = nrow(nn_idx))
      diff_mat    <- nn_clusters != cluster_labels
      frac_diff   <- rowMeans(diff_mat)
    }

    # periphery_score: outer-edge of the full point cloud (all neighbours)
    periphery_score <- numeric(n_cells)

    if (!is.null(coords) && n_dims == 2L) {
      # 2D angular-gap method
      twopi <- 2 * pi
      for (i in seq_len(n_cells)) {
        nb_idx <- nn_idx[i, ]
        dx <- coords[nb_idx, 1L] - coords[i, 1L]
        dy <- coords[nb_idx, 2L] - coords[i, 2L]
        ang <- sort(atan2(dy, dx))
        gaps <- diff(c(ang, ang[1L] + twopi))
        periphery_score[i] <- max(gaps) / twopi
      }
    } else if (!is.null(coords)) {
      # 3D resultant
      for (i in seq_len(n_cells)) {
        nb_idx <- nn_idx[i, ]
        d  <- sweep(coords[nb_idx, , drop = FALSE], 2L, coords[i, ])
        nrm <- sqrt(rowSums(d^2))
        nrm[nrm == 0] <- 1
        u  <- d / nrm
        periphery_score[i] <- sqrt(sum(colSums(u)^2)) / length(nb_idx)
      }
    } else {
      periphery_score <- rep(NA_real_, n_cells)
    }

    border_score <- pmax(frac_diff, periphery_score, na.rm = TRUE)
    if (all(is.na(periphery_score))) border_score <- frac_diff

    border_results[[slot]] <- data.frame(
      frac_diff       = frac_diff,
      periphery_score = periphery_score,
      border_score    = border_score,
      nn_mean_dist    = nb$nn_mean_dist,
      cluster         = cluster_labels,
      stringsAsFactors = FALSE
    )
  }

  # -- High-dimensional reference (frac_diff only) ---------------------------
  hd <- neighbors$high_dim
  if (!is.null(hd)) {
    if (!is.null(hd$nn_idx)) {
      kh <- ncol(hd$nn_idx)
      nn_clusters_hd <- matrix(cluster_labels[hd$nn_idx], nrow = nrow(hd$nn_idx))
      frac_diff_hd <- rowMeans(nn_clusters_hd != cluster_labels)
      deg <- rep(kh, n_cells)
    } else {
      frac_diff_hd <- numeric(n_cells)
      deg <- tabulate(hd$nn_i, nbins = n_cells)
      inter <- tabulate(hd$nn_i[cluster_labels[hd$nn_j] != cluster_labels[hd$nn_i]],
                        nbins = n_cells)
      frac_diff_hd <- ifelse(deg > 0, inter / deg, 0)
    }
    border_results[["high_dim"]] <- data.frame(
      frac_diff       = frac_diff_hd,
      periphery_score = NA_real_,
      border_score    = frac_diff_hd,
      nn_mean_dist    = hd$nn_mean_dist,
      cluster         = cluster_labels,
      stringsAsFactors = FALSE
    )
  }

  fcs_join_obj$reduction_borders <- border_results

  if (!"object_history" %in% names(fcs_join_obj)) {
    message("Consider running FCSimple::fcs_audit() on the object.")
  }
  try(expr = fcs_join_obj[["object_history"]] <- append(
    fcs_join_obj[["object_history"]],
    paste0("reduction_borders on ", tolower(algorithm), " (",
           paste(names(border_results), collapse = ", "), "): ", Sys.time())
  ), silent = TRUE)

  message("Border scores computed for ", length(border_results), " slot(s).")
  return(fcs_join_obj)
}