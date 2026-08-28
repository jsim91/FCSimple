# utils_colors.R
# Shared color palette and cluster color optimization utilities used by
# fcs_plot_reduction() and fcs_plot_reduction_3d().

# -- Perceptually-balanced categorical palette -----------------------------------
# Uses RColorBrewer Set2 (<=8), Paired (<=12), then farthest-point sampling
# in CIE L*a*b* over an HCL candidate grid for larger n.
fcs_color_palette <- function(n) {
  if (n <= 8L)  return(RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)])
  if (n <= 12L) return(RColorBrewer::brewer.pal(max(3L, n), "Paired")[seq_len(n)])

  # Farthest-point sampling in CIE Lab over an HCL candidate grid
  h_vals <- seq(0, 355, by = 5)
  c_vals <- seq(40, 100, by = 10)
  l_vals <- seq(40, 85,  by = 5)
  grid   <- expand.grid(H = h_vals, C = c_vals, L = l_vals)
  hex    <- grDevices::hcl(h = grid$H, c = grid$C, l = grid$L, fixup = FALSE)
  ok     <- !is.na(hex)
  hex    <- hex[ok]
  if (length(hex) < n) return(scales::hue_pal()(n))

  rgb_m <- t(grDevices::col2rgb(hex)) / 255
  lab   <- grDevices::convertColor(rgb_m, from = "sRGB", to = "Lab")

  sel      <- integer(n)
  sel[1L]  <- which.max(lab[, 2L]^2 + lab[, 3L]^2)   # seed: most chromatic
  min_d    <- sqrt(rowSums(sweep(lab, 2, lab[sel[1L], ])^2))
  for (i in seq.int(2L, n)) {
    sel[i] <- which.max(min_d)
    min_d   <- pmin(min_d, sqrt(rowSums(sweep(lab, 2, lab[sel[i], ])^2)))
  }
  hex[sel]
}

# -- Space-aware cluster color assignment (DSatur) -------------------------------
# Builds a k-NN graph on cluster centroids (k=4, or fewer for small n) using
# Euclidean distance, then assigns colors from base_colors via the DSatur
# (degree-of-saturation) greedy algorithm: at each step, pick the uncolored
# cluster with the most already-colored neighbors (breaking ties by minimum
# spatial distance to any colored cluster), and assign the unused palette color
# whose minimum CIE L*a*b* perceptual distance to already-colored neighbors is
# greatest.  This maximizes local contrast so that bordering clusters are easily
# distinguishable.
#
# @param centroids A numeric matrix of cluster centroids:
#   rows = clusters, cols = dimensions (2 for 2D, 3 for 3D, etc.)
#   rownames are used as cluster names.
# @param base_colors Character vector of hex colors (one per cluster).
# @param k Integer; k for the k-NN graph (default 4).
# @return Named character vector (names = cluster names, values = hex colors).
# -- Reduction slot resolver -----------------------------------------------------
# Given a reduction name (e.g. "UMAP" or "umap") and required dimensionality,
# resolves to the correct slot name in fcs_join_obj.  Looks for dimensionalized
# names first (umap_2d, umap_3d), falls back to legacy flat names (umap).
# @param obj fcs_join_obj list
# @param reduction character; case-insensitive reduction name
# @param dims integer; 2 or 3 — exact dimensionality required
# @return character; the actual slot name in obj
.resolve_reduction_slot <- function(obj, reduction, dims) {
  red <- tolower(reduction)
  target <- paste0(red, "_", dims, "d")
  if (target %in% names(obj)) return(target)
  # Legacy fallback for objects created before dimensional slots
  if (red %in% names(obj)) {
    coords <- obj[[red]][["coordinates"]]
    if (!is.null(coords) && ncol(coords) >= dims) return(red)
  }
  available <- intersect(
    c("umap_2d", "umap_3d", "tsne_2d", "tsne_3d", "umap", "tsne"),
    names(obj)
  )
  stop("No ", dims, "D reduction for '", reduction, "' found. Available: ",
       paste(available, collapse = ", "))
}

# -- Primary reduction & shared cluster colour mapping ---------------------------
# Returns the first-created reduction slot name.  This drives cluster colouring
# so that 2D and 3D plots of the same clusters use identical colour assignments
# regardless of which embedding is being drawn.
.get_primary_reduction_slot <- function(obj) {
  if (!is.null(obj$reduction_order) && length(obj$reduction_order) > 0L)
    return(obj$reduction_order[1L])
  reds <- grep("^(umap|tsne)(_[23]d)?$", names(obj), value = TRUE)
  if (length(reds) > 0L) return(reds[1L])
  stop("No reduction embedding found in the object.")
}

# Computes a deterministic cluster -> colour mapping using DSatur on a chosen
# reduction's cluster centroids.  Colouring from the primary reduction ensures
# consistent colours across 2D and 3D plots.
# @param obj fcs_join_obj list
# @param cluster_numbers character vector (1 per cell); may be substituted labels
# @param reduction_slot character; which reduction coordinates to build centroids from
# @return named character vector (cluster label -> hex colour)
.compute_cluster_color_mapping <- function(obj, cluster_numbers, reduction_slot) {
  coords <- obj[[reduction_slot]]$coordinates
  uclus  <- unique(cluster_numbers)
  uclus  <- uclus[order(uclus)]
  centroids <- matrix(NA_real_, nrow = length(uclus), ncol = ncol(coords))
  rownames(centroids) <- uclus
  for (i in seq_along(uclus)) {
    idx <- which(cluster_numbers == uclus[i])
    centroids[i, ] <- apply(coords[idx, , drop = FALSE], 2L, stats::median)
  }
  base_colors <- fcs_color_palette(length(uclus))
  optimize_cluster_colors(centroids, base_colors)
}

optimize_cluster_colors <- function(centroids, base_colors, k = 4L) {
  cluster_names <- rownames(centroids)
  n_clusters    <- nrow(centroids)
  n_colors      <- length(base_colors)

  # Euclidean distance matrix between cluster centroids
  dist_matrix <- as.matrix(stats::dist(centroids, method = "euclidean"))

  # Build k-NN graph (k = min(k, n_clusters - 1))
  k_eff <- min(k, n_clusters - 1L)
  neighbors <- t(apply(dist_matrix, 1L, function(d) order(d)[seq_len(k_eff + 1L)][-1L]))
  rownames(neighbors) <- cluster_names

  # Convert base_colors to CIE Lab for perceptual distance calculations
  rgb_m    <- t(grDevices::col2rgb(base_colors)) / 255
  lab_cols <- grDevices::convertColor(rgb_m, from = "sRGB", to = "Lab")
  rownames(lab_cols) <- seq_len(n_colors)

  color_assignment <- integer(n_clusters)
  names(color_assignment) <- cluster_names
  color_used <- logical(n_colors)

  # Seed: cluster closest to the global centroid
  centroid_global <- colMeans(centroids)
  start_cluster <- names(which.min(rowSums(sweep(centroids, 2, centroid_global)^2)))
  color_assignment[start_cluster] <- 1L
  color_used[1L] <- TRUE
  colored_set <- start_cluster

  while (length(colored_set) < n_clusters) {
    uncolored <- setdiff(cluster_names, colored_set)

    # DSatur: count already-colored neighbors for each uncolored cluster
    saturation <- sapply(uncolored, function(uc) {
      sum(neighbors[uc, ] %in% colored_set)
    })

    # Break ties by minimum distance to any already-colored cluster
    min_dist_to_colored <- sapply(uncolored, function(uc) {
      min(dist_matrix[uc, colored_set])
    })

    # Pick the uncolored cluster with most colored neighbors;
    # ties broken by smallest minimum distance to a colored cluster
    best_idx <- which.max(saturation - min_dist_to_colored / max(dist_matrix))
    next_cluster <- uncolored[best_idx]

    # Gather colors already used by this cluster's already-colored neighbors
    nbr_colored <- intersect(neighbors[next_cluster, ], colored_set)
    nbr_col_idx <- color_assignment[nbr_colored]

    available <- which(!color_used)

    if (length(available) > 0L && length(nbr_col_idx) > 0L) {
      # Choose the available color whose minimum Lab distance
      # to any neighbor color is maximized
      scores <- sapply(available, function(ci) {
        min(sapply(nbr_col_idx, function(ni) {
          sqrt(sum((lab_cols[ci, ] - lab_cols[ni, ])^2))
        }))
      })
      best_color <- available[which.max(scores)]
    } else if (length(available) > 0L) {
      best_color <- available[1L]
    } else {
      best_color <- 1L  # all colors used (edge case)
    }

    color_assignment[next_cluster] <- best_color
    color_used[best_color] <- TRUE
    colored_set <- c(colored_set, next_cluster)
  }

  result        <- base_colors[color_assignment]
  names(result) <- cluster_names
  result
}