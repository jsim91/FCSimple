#' @title Compute Local Neighbor Graphs and r-Hop Neighborhoods
#'
#' @description
#'   Runs a k-nearest-neighbor search on every dimensionality reduction
#'   embedding present in the object (UMAP/t-SNE, 2D/3D) and resolves the
#'   high-dimensional reference neighbor graph from clustering.  For each
#'   space it additionally computes \emph{r-hop} neighborhoods by graph
#'   traversal over the k-NN graph, so that a cell's semi-local neighborhood
#'   can be grown beyond its immediate k neighbours without ever re-querying
#'   raw coordinates.
#'
#' @param fcs_join_obj
#'   A list produced by \code{FCSimple::fcs_reduce_dimensions()} and, for the
#'   high-dimensional reference, \code{FCSimple::fcs_cluster()}.
#'
#' @param k
#'   Integer; number of nearest neighbours.  Default \code{30}.
#'
#' @param n_hops
#'   Integer >= 1; radius, in graph hops, of the semi-local neighborhood
#'   stored for each cell.  A value of 1 stores only immediate k-NN
#'   neighbours; larger values expand neighborhoods by traversing the k-NN
#'   graph (reusing the high-dimensional adjacency matrix when available).
#'   Default \code{1}.
#'
#' @param high_k
#'   Integer or \code{NULL}.  When \code{NULL} (default), the high-dimensional
#'   reference reuses the clustering graph stored on the object (the
#'   \code{adjacency_matrix} or \code{search} result from
#'   \code{fcs_cluster()}).  When set, the high-dimensional k-NN graph is
#'   instead \emph{re-derived} from the clustering data matrix at this many
#'   neighbours, guaranteeing parity with \code{k} when the two are matched.
#'
#' @param high_use_rep
#'   Character; which high-dimensional representation to use when
#'   \code{high_k} is set.  One of \code{"data"} (default; uses
#'   batch-corrected data if present, else raw \code{data}),
#'   \code{"pca"}.  Ignored when \code{high_k} is \code{NULL}.
#'
#' @param num_cores
#'   Integer; number of CPU cores for the parallel neighbour search.
#'   Default \code{ceiling(parallel::detectCores() / 2)}.
#'
#' @details
#'   \strong{Embeddings:} \code{"umap_2d"}, \code{"umap_3d"},
#'   \code{"tsne_2d"}, and \code{"tsne_3d"} (plus legacy \code{"umap"} /
#'   \code{"tsne"}) are auto-discovered.  For each, \code{RANN::nn2()} is run
#'   using a multisession \code{future} backend (with sequential fallback) and
#'   the k nearest-neighbour indices and distances (excluding self) are stored.
#'
#'   \strong{High-dimensional reference:} the existing clustering graph is
#'   reused.  If \code{fcs_join_obj$search} is present its neighbour columns
#'   are used; otherwise \code{fcs_join_obj$adjacency_matrix} is summarised via
#'   \code{Matrix::summary()}.
#'
#'   \strong{r-hop neighborhoods:} for every space, an undirected k-NN graph is
#'   built and \code{igraph::neighborhood(order = n_hops, mindist = 1)} returns
#'   the set of cells reachable within \code{n_hops}.  These sets drive
#'   \code{fcs_neighborhood_consistency()} and the semi-local border scoring,
#'   and are matched in definition and cardinality across high and low
#'   dimensional spaces.
#'
#' @return
#'   The input \code{fcs_join_obj} augmented with \code{$neighbors}, a named
#'   list with one entry per embedding (e.g. \code{"umap_2d"}) plus a
#'   \code{"high_dim"} entry.  Each embedding entry contains:
#'   \describe{
#'     \item{\code{nn_idx}}{Numeric n x k matrix of neighbour indices.}
#'     \item{\code{nn_dist}}{Numeric n x k matrix of neighbour distances.}
#'     \item{\code{nn_mean_dist}}{Numeric vector of per-cell mean neighbour
#'       distances.}
#'     \item{\code{hop_idx}}{List of length n; each element is an integer
#'       vector of r-hop neighbour cell indices.}
#'   }
#'   The \code{"high_dim"} entry additionally stores edge-list columns
#'   \code{nn_i}/\code{nn_j}, an r-hop neighbour list \code{hop_idx}, and the
#'   \code{features} the high-dimensional graph was derived from, or is
#'   \code{NULL} when no clustering graph is available.
#'
#'   \code{object_history} is appended with a timestamped entry.
#'
#' @examples
#' \dontrun{
#'   obj <- FCSimple::fcs_reduce_dimensions(obj, algorithm = "umap")
#'   obj <- FCSimple::fcs_embedding_neighbors(obj, k = 30, n_hops = 2)
#'   length(obj$neighbors$umap_2d$hop_idx[[1]])
#' }
#'
#' @seealso
#'   FCSimple::fcs_neighborhood_consistency,
#'   FCSimple::fcs_reduction_borders
#'
#' @importFrom RANN nn2
#' @importFrom parallel detectCores
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom Matrix summary
#' @importFrom igraph make_empty_graph add_edges neighborhood
#' @export
fcs_embedding_neighbors <- function(
    fcs_join_obj,
    k = 30,
    n_hops = 1,
    high_k = NULL,
    high_use_rep = "data",
    num_cores = ceiling(parallel::detectCores() / 2))
{
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("Package 'RANN' is required. Install with: install.packages('RANN')")
  if (!requireNamespace("parallel", quietly = TRUE))
    stop("Package 'parallel' is required.")
  if (!requireNamespace("future", quietly = TRUE))
    stop("Package 'future' is required. Install with: install.packages('future')")
  if (!requireNamespace("future.apply", quietly = TRUE))
    stop("Package 'future.apply' is required. Install with: install.packages('future.apply')")
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required. Install with: install.packages('igraph')")

  n_hops <- as.integer(n_hops)
  if (n_hops < 1L) stop("'n_hops' must be >= 1.")

  # -- Discover available embeddings -----------------------------------------
  emb_candidates <- grep("^(umap|tsne)(_[23]d)?$", names(fcs_join_obj),
                         value = TRUE)
  if (length(emb_candidates) == 0L)
    stop("No reduction embeddings found in fcs_join_obj.  ",
         "Run fcs_reduce_dimensions() first.")

  embedding_list <- list()
  for (red in emb_candidates) {
    coords <- fcs_join_obj[[red]]$coordinates
    if (is.null(coords)) next
    nd <- ncol(coords)
    if (!nd %in% c(2L, 3L)) {
      warning("Skipping '", red, "': ", nd,
              " dimensions (only 2D and 3D are supported).")
      next
    }
    slot_name <- if (grepl("_\\d+d$", red)) red else paste0(red, "_", nd, "d")
    embedding_list[[slot_name]] <- coords
  }
  n_cells <- if (length(embedding_list) > 0L) nrow(embedding_list[[1]]) else NULL
  if (is.null(n_cells))
    stop("No valid 2D/3D embeddings found.")

  # -- Resolve core count ----------------------------------------------------
  n_cores <- num_cores
  if (n_cores > parallel::detectCores()) {
    warning(n_cores, " cores specified but only ", parallel::detectCores(),
            " available.  Proceeding with max available cores.")
    n_cores <- parallel::detectCores()
  } else if (n_cores == 0L) {
    n_cores <- parallel::detectCores()
  }
  if (n_cores > 16L) {
    warning("Using ", n_cores,
            " cores may cause excessive overhead.  Consider 8-16 cores for optimal performance.")
  }

  options(future.globals.maxSize = Inf)
  future::plan("multisession", workers = n_cores)
  on.exit(future::plan(future::sequential), add = TRUE)

  num_neighbors <- k + 1L  # +1 because self will be excluded

  # -- Helper: build undirected graph and r-hop neighbour lists --------------
  build_hop_sets <- function(ii, jj, n, hops) {
    g <- igraph::make_empty_graph(n = n, directed = FALSE)
    if (length(ii) > 0L) {
      g <- igraph::add_edges(g, edges = as.vector(rbind(ii, jj)))
    }
    hop <- igraph::neighborhood(g, order = hops, mindist = 1L)
    lapply(hop, as.integer)
  }

  neighbors <- list()

  # -- Process each embedding ------------------------------------------------
  for (embed_name in names(embedding_list)) {
    message("Processing ", embed_name, " (", ncol(embedding_list[[embed_name]]),
            "D, k=", k, ", n_hops=", n_hops, ")...")
    coords <- as.matrix(embedding_list[[embed_name]])

    # Split coordinates into chunks for parallel query
    split_sums <- round(seq(from = 1, to = nrow(coords),
                            length.out = n_cores + 1), 0)
    sub_data <- vector("list", n_cores)
    for (i in seq_along(sub_data)) {
      if (i == 1L) {
        sub_data[[i]] <- coords[1:(split_sums[i + 1]), , drop = FALSE]
      } else {
        sub_data[[i]] <- coords[(split_sums[i] + 1):(split_sums[i + 1]), ,
                                drop = FALSE]
      }
    }

    # Parallel nn2 with fallback
    search_out <- tryCatch({
      future.apply::future_lapply(sub_data, FUN = function(x) {
        RANN::nn2(data = coords, query = x, k = num_neighbors,
                  treetype = "kd", searchtype = "standard")
      })
    }, error = function(e) {
      warning("multisession future_lapply failed: ", conditionMessage(e),
              ".  Falling back to sequential search.")
      NULL
    })

    if (is.null(search_out)) {
      future::plan(future::sequential)
      search_out <- lapply(sub_data, FUN = function(x) {
        RANN::nn2(data = coords, query = x, k = num_neighbors,
                  treetype = "kd", searchtype = "standard")
      })
    }

    nn_idx  <- do.call(rbind, lapply(search_out, `[[`, 1L))
    nn_idx  <- nn_idx[, 2:num_neighbors, drop = FALSE]
    nn_dist <- do.call(rbind, lapply(search_out, `[[`, 2L))
    nn_dist <- nn_dist[, 2:num_neighbors, drop = FALSE]

    storage.mode(nn_idx) <- "integer"

    ii <- rep(seq_len(n_cells), each = k)
    jj <- as.vector(t(nn_idx))

    neighbors[[embed_name]] <- list(
      nn_idx       = nn_idx,
      nn_dist      = nn_dist,
      nn_mean_dist = rowMeans(nn_dist),
      hop_idx      = build_hop_sets(ii, jj, n_cells, n_hops)
    )
  }

  # -- Resolve high-dimensional reference ------------------------------------
  # By default reuse the clustering graph; optionally re-derive at `high_k`
  # from the clustering data matrix for explicit high/Low-D cardinality parity.
  high_dim <- NULL
  hd_i <- integer(0L)
  hd_j <- integer(0L)

  resolve_high_data <- function(obj, use_rep) {
    if ("batch_correction" %in% names(obj) &&
        !is.null(obj$batch_correction$data)) {
      message("Using batch-corrected data as the high-dimensional reference.")
      return(as.matrix(obj$batch_correction$data))
    }
    if (identical(use_rep, "pca") && "pca" %in% names(obj) &&
        !is.null(obj$pca$pca_data)) {
      message("Using PCA data as the high-dimensional reference.")
      return(as.matrix(obj$pca$pca_data))
    }
    if ("data" %in% names(obj) && !is.null(obj$data)) {
      return(as.matrix(obj$data))
    }
    stop("No high-dimensional data matrix available to re-derive the ",
         "reference graph by 'high_k'. Provide 'data', 'pca$pca_data', or ",
         "'batch_correction$data'.")
  }

  if (!is.null(high_k)) {
    high_k <- as.integer(high_k)
    if (high_k < 1L) stop("'high_k' must be >= 1.")
    high_data <- resolve_high_data(fcs_join_obj, high_use_rep)
    if (nrow(high_data) != n_cells)
      stop("High-dimensional reference has ", nrow(high_data),
           " rows but embeddings have ", n_cells, " cells.")
    kh <- min(high_k, n_cells - 1L)
    message("Re-deriving high-dimensional reference at k=", kh, ".")
    hd_nn <- RANN::nn2(data = high_data, query = high_data, k = kh + 1L)
    nn_idx_high <- hd_nn$nn.idx[, 2:(kh + 1L), drop = FALSE]
    nn_dist_high <- hd_nn$nn.dists[, 2:(kh + 1L), drop = FALSE]
    storage.mode(nn_idx_high) <- "integer"
    hd_i <- rep(seq_len(n_cells), each = kh)
    hd_j <- as.vector(t(nn_idx_high))
    high_dim <- list(
      features     = .fcs_features(high_data),
      nn_idx       = nn_idx_high,
      nn_i         = hd_i,
      nn_j         = hd_j,
      nn_dist      = nn_dist_high,
      nn_mean_dist = rowMeans(nn_dist_high),
      hop_idx      = build_hop_sets(hd_i, hd_j, n_cells, n_hops)
    )
  } else if (!is.null(fcs_join_obj$search)) {
    nn_idx_high <- fcs_join_obj$search[
      , 2:min(ncol(fcs_join_obj$search), k + 1L), drop = FALSE]
    kh <- ncol(nn_idx_high)
    storage.mode(nn_idx_high) <- "integer"
    hd_i <- rep(seq_len(n_cells), each = kh)
    hd_j <- as.vector(t(nn_idx_high))
    high_dim <- list(
      features     = fcs_join_obj$search_features,
      nn_idx       = nn_idx_high,
      nn_i         = hd_i,
      nn_j         = hd_j,
      nn_dist      = NULL,
      nn_mean_dist = rep(NA_real_, n_cells),
      hop_idx      = build_hop_sets(hd_i, hd_j, n_cells, n_hops)
    )
  } else if (!is.null(fcs_join_obj$adjacency_matrix)) {
    coords <- Matrix::summary(fcs_join_obj$adjacency_matrix)
    hd_i <- coords$i
    hd_j <- coords$j
    high_dim <- list(
      features     = fcs_join_obj$search_features,
      nn_idx       = NULL,
      nn_i         = hd_i,
      nn_j         = hd_j,
      nn_dist      = NULL,
      nn_mean_dist = rep(NA_real_, n_cells),
      hop_idx      = build_hop_sets(hd_i, hd_j, n_cells, n_hops)
    )
  }
  neighbors[["high_dim"]] <- high_dim

  fcs_join_obj$neighbors <- neighbors

  if (!"object_history" %in% names(fcs_join_obj)) {
    message("Consider running FCSimple::fcs_audit() on the object.")
  }
  try(expr = fcs_join_obj[["object_history"]] <- append(
    fcs_join_obj[["object_history"]],
    paste0("embedding_neighbors (", paste(names(embedding_list), collapse = ", "),
           if (is.null(high_dim)) "" else ", high_dim", "; n_hops=", n_hops,
           "): ", Sys.time())
  ), silent = TRUE)

  message("Neighbour graphs and r-hop sets computed for ",
          length(embedding_list), " embedding(s) and ",
          if (is.null(high_dim)) "no" else "a",
          " high-dimensional reference.")
  return(fcs_join_obj)
}