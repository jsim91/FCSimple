#’ @title Initialize Gating Object
#’
#’ @description
#’   Wraps a joined FCSimple object into a gating‐friendly structure.  
#’   Records the raw data matrix, sample sources, optional run dates,  
#’   and initializes an empty gating history.  
#’
#’ @param fcs_obj  
#’   A list‐like FCSimple object (e.g. from FCSimple::fcs_join) containing at  
#’   least:  
#’   - data: numeric matrix or data.frame of events × channels  
#’   - source: character vector of sample identifiers  
#’   - run_date: (optional) acquisition dates  
#’
#’ @return  
#’   An object of class “fcs_gating_object” with elements:  
#’   - data  
#’   - source  
#’   - run_date (if present)  
#’   - object_history: character vector recording the gating steps  
#’
#’ @seealso  
#’   FCSimple::fcs_join, FCSimple::fcs_update  
#’
#’ @importFrom FCSimple fcs_update  
#’ @export
fcs_gating_object <- function(fcs_obj) {
  print(paste0("Initiating an object for gating with the assumption that data to be gated was collected using fluorescence cytometry. Cytof is not supported (currently)."))
  print("Data for gating is assumed to be in slot fcs_obj[['data']].")
  list_obj <- list(data = fcs_obj$data, 
                   source = fcs_obj$source)
  if('run_date' %in% names(fcs_obj)) {
    list_obj[['run_date']] <- fcs_obj$run_date
  }
  list_obj <- FCSimple::fcs_update(fcs_join_obj = list_obj, instrument_type = 'flow')
  class(list_obj) <- 'fcs_gating_object'
  return(list_obj)
}

#’ @title Gate Cells by SSC-A vs FSC-A
#’
#’ @description
#’   Identifies the main cell population by fitting a bivariate Gaussian mixture  
#’   model (Mclust) on side‐scatter (SSC‐A) vs forward‐scatter (FSC‐A), drops the  
#’   nearest/farthest mixture components, computes a convex hull of retained events,  
#’   and returns the gating mask (and adds it to the gating object if provided).  
#’
#’ @param object  
#’   An “fcs_gating_object”, or a data.frame/matrix with columns `side_a` and  
#’   `forward_a`.  
#’ @param tree_name  
#’   Character; name of the gating tree (default “tree1”).  
#’ @param parent_name  
#’   Character; name of the parent gate in the tree (default “none”).  
#’ @param gate_name  
#’   Character; new gate label (no “+”/“–” allowed).  
#’ @param side_a  
#’   Character; name of side‐scatter channel (default “SSC-A”).  
#’ @param forward_a  
#’   Character; name of forward‐scatter channel (default “FSC-A”).  
#’ @param drop_nearest  
#’   Integer; number of nearest‐to‐origin components to drop (default 1).  
#’ @param drop_farthest  
#’   Integer; number of farthest‐from‐origin components to drop (default 1).  
#’
#’ @details
#’ - Fits `mclust::Mclust(df[,c(side_a, forward_a)], modelNames="VVV")`.  
#’ - Computes center distances to origin, orders components, drops extremes.  
#’ - Builds convex hull of retained points with `chull()`.  
#’ - Binary mask: event inside hull = 1, else 0.  
#’
#’ @return  
#’ If `object` is an “fcs_gating_object”, returns it with:  
#’ - `object$gate_trees[[tree_name]][[gate_name]]`: list with  
#’   - mask, tree, parent, feature_side_a, feature_forward_a, gate_fn,  
#’     chull_selection (hull vertices + cluster indices).  
#’ Otherwise returns a list with:  
#’ - cells: data.frame of gated events  
#’ - hull_vertices: data.frame of hull points  
#’ - kept_clusters: integer indices  
#’ - dropped_clusters: integer indices  
#’
#’ @examples
#’ \dontrun{
#’   go <- fcs_gating_object(joined)
#’   gated <- fcs_gate_cells(go, gate_name = "cells")
#’ }
#’
#’ @seealso  
#’   fcs_gating_object, mclust::Mclust, sp::point.in.polygon  
#’
#’ @importFrom mclust Mclust  
#’ @importFrom sp point.in.polygon  
#’ @export
fcs_gate_cells <- function(object,
                           tree_name = 'tree1', 
                           parent_name = 'none', 
                           gate_name = 'cells', 
                           side_a = "SSC-A",
                           forward_a = "FSC-A",
                           # downsample_and_infer_gate_members = FALSE, # implement method (and benchmark run time) that allows for RANN knn=1? to identify gate membership from downsampled data; there should be an existing package for this.
                           drop_nearest = 1,
                           drop_farthest = 1) {
  require(mclust)
  require(sp)
  if(regexpr(pattern = '(\\+|\\-)', text = gate_name)!=-1) {
    stop("'gate_name' should not include '+' or '-'")
  }
  if(class(object)=='fcs_gating_object') {
    print("General note: if this gating tree borrows from another, consider copying the common gate path to this new tree before continuing. Each tree is linear and is independent of all other trees.")
    if(!'gate_trees' %in% names(object)) {
      object[['gate_trees']] <- list()
    }
    if(!tree_name %in% names(object[['gate_trees']])) {
      object[['gate_trees']][[tree_name]] <- list()
    }
    tmpdf <- as.data.frame(object[['data']]); row.names(tmpdf) <- 1:nrow(tmpdf)
    if(parent_name!='none') {
      warning("It's recommended that when using 'fcs_gate_cells' that you start from the root data--no preceeding gates (the 'parent_name' default). Proceed with caution.")
      if(!parent_name %in% object[['gate_trees']][[tree_name]]) {
        stop("Could not find 'parent_name' in specified 'tree_name'. If 'parent_name' is not 'none', the specified 'parent_name' must exist in the specified 'tree_name'.")
      } else {
        parent_index <- which(names(object[['gate_trees']][[tree_name]])==parent_name)
        mask_list <- lapply(X = object[['gate_trees']][[tree_name]][1:parent_index], FUN = function(arg) return(arg[['mask']]))
        cell_mask <- do.call(pmin, mask_list)
        df <- tmpdf[cell_mask==1,] # 1 = in parent gate; 0 = not in parent gate
        names(cell_mask) <- row.names(tmpdf)
      }
    } else {
      df <- tmpdf; rm(tmpdf)
    }
  } else if(class(object)=='data.frame') {
    df <- object
  } else if('matrix' %in% class(object)){
    df <- as.data.frame(object)
  } else {
    stop("object should be of class: data.frame, matrix, or fcs_gating_object")
  }
  if(!all(c(side_a, forward_a) %in% names(df))) {
    stop("df must contain columns '", side_a, "' and '", forward_a, "'")
  }
  
  mcl <- Mclust(df[,c(side_a, forward_a)], modelNames = "VVV")
  
  ctrs <- as.data.frame(t(mcl$parameters$mean))
  names(ctrs) <- c(side_a, forward_a)
  ctrs$dist <- sqrt(ctrs[[side_a]]^2 + ctrs[[forward_a]]^2)
  
  G <- nrow(ctrs)
  
  if(drop_nearest + drop_farthest >= G) {
    stop("drop_nearest + drop_farthest must be less than the total number of components (", G, ")")
  }
  if(drop_nearest < 0 || drop_farthest < 0) {
    stop("drop_nearest and drop_farthest must be non-negative integers")
  }
  
  ord <- order(ctrs$dist)
  
  to_drop <- integer(0)
  if(drop_nearest  > 0) to_drop <- c(to_drop, ord[seq_len(drop_nearest)])
  if(drop_farthest > 0) to_drop <- c(to_drop, ord[(G - drop_farthest + 1L):G])
  
  keep_clusters <- setdiff(seq_len(G), to_drop)
  
  cls <- mcl$classification
  good_idx <- which(cls %in% keep_clusters)
  good_df  <- df[good_idx, , drop=FALSE]
  
  hull_i <- chull(good_df[[side_a]], good_df[[forward_a]])
  hull_pts <- good_df[hull_i, c(side_a, forward_a)]
  
  pip <- point.in.polygon(df[[side_a]], 
                          df[[forward_a]],
                          hull_pts[[side_a]], 
                          hull_pts[[forward_a]])
  inside <- pip > 0
  gated  <- df[inside, , drop=FALSE]
  
  if(parent_name!='none') {
    current_gate <- as.numeric(inside); names(current_gate) <- row.names(df)
    pos <- match(names(current_gate), names(cell_mask))
    final_mask <- cell_mask
    final_mask[pos] <- current_gate
    names(final_mask) <- NULL
  } else {
    final_mask <- as.numeric(inside)
  }
  if(class(object)=='fcs_gating_object') {
    new_branch <- list('mask' = final_mask, # T,F to 1,0
                       'tree' = tree_name, 
                       'parent' = parent_name, 
                       'feature_side_a' = side_a, 
                       'feature_forward_a' = forward_a, 
                       'gate_fn' = 'fcs_gate_cells', 
                       'chull_selection' = list('hull_Vertices' = hull_pts, 
                                                'kept_clusters' = keep_clusters, 
                                                'dropped_clusters' = to_drop))
    object[['gate_trees']][[tree_name]] <- append(object[['gate_trees']][[tree_name]], list(new_branch))
    names(object[['gate_trees']][[tree_name]])[length(object[['gate_trees']][[tree_name]])] <- gate_name
    return(object)
  } else {
    outlist <- list(cells = gated,
                    hull_vertices = hull_pts,
                    kept_clusters = keep_clusters,
                    dropped_clusters = to_drop)
    return(outlist)
  }
}

#’ @title Plot Gated Cells on SSC vs FSC
#’
#’ @description
#’   Scatter‐plots all events in gray and gated events in blue, overlaying  
#’   the convex hull in red.  
#’
#’ @param object  
#’   An “fcs_gating_object” containing `gate_trees[[gate_tree]][[gate_name]]`.  
#’ @param gate_tree  
#’   Character; name of the gating tree (default “tree1”).  
#’ @param gate_name  
#’   Character; name of the gate branch to plot (default “cells”).  
#’ @param downsample_size  
#’   Integer; maximum points to plot (default 100000).  
#’
#’ @return  
#’   A ggplot2 object.  
#’
#’ @examples
#’ \dontrun{
#’   go <- fcs_gate_cells(gating_obj, gate_name = "cells")
#’   p <- fcs_plot_cells(go, gate_name = "cells")
#’   print(p)
#’ }
#’
#’ @seealso  
#’   fcs_gate_cells, ggplot2::ggplot  
#’
#’ @importFrom ggplot2 ggplot aes geom_point geom_polygon theme_bw labs theme element_text  
#’ @importFrom rlang sym  
#’ @export
fcs_plot_cells <- function(object, 
                           gate_tree, 
                           gate_name = 'cells', 
                           downsample_size = 100000) {
  require(ggplot2)
  require(rlang)
  
  if(class(object)!='fcs_gating_object') {
    stop("input object should be of class: fcs_gating_object")
  }
  gate_res = fcs_gate_obj$gate_trees[[gate_tree]][[gate_name]][['chull_selection']]
  side_a <- fcs_gate_obj$gate_trees[[gate_tree]][[gate_name]]$feature_side_a
  forward_a <- fcs_gate_obj$gate_trees[[gate_tree]][[gate_name]]$feature_forward_a
  cells_all <- fcs_gate_obj$data[,c(side_a, forward_a)]
  cells <- fcs_gate_obj$data[,c(side_a, forward_a)][fcs_gate_obj$gate_trees[[gate_tree]][[gate_name]]$mask==1,]
  hull <- fcs_gate_obj$gate_trees[[gate_tree]][[gate_name]][['chull_selection']]$hull_Vertices
  
  x_sym <- sym(side_a)
  y_sym <- sym(forward_a)
  
  if(!is.na(downsample_size)) {
    ds_scaled <- ceiling((downsample_size/nrow(cells_all))*nrow(cells))
    if(nrow(cells_all)>downsample_size) {
      cells_all <- cells_all[sample(1:nrow(cells_all),downsample_size,replace=F),]
    }
    if(nrow(cells)>downsample_size) {
      cells <- cells[sample(1:nrow(cells),ds_scaled,replace=F),]
    }
  }
  
  pl <- ggplot(cells_all, aes(x = !!x_sym, y = !!y_sym)) + 
    geom_point(color = "grey20", size = 0.8, alpha = 0.5) + 
    geom_point(data = cells,
               aes(x = !!x_sym, y = !!y_sym),
               color = "blue", size = 1, alpha = 0.5) +
    geom_polygon(data = hull,
                 aes(x = !!x_sym, y = !!y_sym),
                 fill  = NA,
                 color = "red", linewidth = 1.5) + 
    theme_bw() +
    labs(x = side_a, y = forward_a) + 
    theme(axis.title = element_text(size = 20, face = 'bold'), 
          axis.text = element_text(size = 14), 
          plot.title = element_text(size = 24, hjust = 0.5, face = 'bold'))
  return(pl)
}

#’ @title Gate Singlets by Height vs Area
#’
#’ @description
#’   Identifies singlet events by fitting a robust or least‐squares line  
#’   to A vs H channels, computing residuals, then applying an EM-based or  
#’   curvature-based cutpoint to select events above the threshold.  
#’
#’ @param object  
#’   An “fcs_gating_object”, or a data.frame/matrix with columns `a` and `h`.  
#’ @param tree_name  
#’   Character; gating tree name (default “tree1”).  
#’ @param parent_name  
#’   Character; parent gate (default “cells”).  
#’ @param gate_name  
#’   Character; new gate label (no “+”/“–” allowed).  
#’ @param a  
#’   Character; area channel name (default “SSC-A”).  
#’ @param h  
#’   Character; height channel name (default “SSC-H”).  
#’ @param cut_method  
#’   Character; “gmm”, “flex”, or “both” (default).  
#’ @param fit_method  
#’   Character; “lm” or “rlm” (default “lm”).  
#’ @param trim_frac  
#’   Numeric; fraction to trim when fitting (default 0.01).  
#’ @param curvature_eps  
#’   Numeric; curvature tolerance for flex-cut (default 0.001).  
#’
#’ @details
#’ - Subsets central data by trim_frac.  
#’ - Fits `lm()` or `MASS::rlm()` of h ~ a.  
#’ - Computes residuals and calls fcs_get_em_cutpoint() for threshold.  
#’ - Masks events above cutpoint as singlets.  
#’
#’ @return  
#’ If `object` is an “fcs_gating_object”, adds branch under  
#’ `gate_trees[[tree_name]][[gate_name]]` with positive/negative masks,  
#’ threshold, fit parameters, and returns the updated object.  
#’ Otherwise returns a data.frame of gated events.  
#’
#’ @examples
#’ \dontrun{
#’   go <- fcs_gating_object(joined)
#’   go2 <- fcs_gate_singlets(go, gate_name = "singlets")
#’ }
#’
#’ @seealso  
#’   fcs_gate_cells, fcs_get_em_cutpoint, MASS::rlm  
#’
#’ @importFrom MASS lm rlm  
#’ @export
fcs_gate_singlets <- function(object,
                              tree_name = 'tree1', 
                              parent_name = 'cells', 
                              gate_name = 'singlets', 
                              a = "SSC-A",
                              h = "SSC-H",
                              cut_method = c('gmm','flex','both'), 
                              fit_method = c("lm","rlm"),
                              trim_frac = 0.01,   # drop extreme 1% when fitting
                              curvature_eps = 0.001
) {
  require(MASS)
  if(regexpr(pattern = '(\\+|\\-)', text = gate_name)!=-1) {
    stop("'gate_name' should not include '+' or '-'")
  }
  if(class(object)=='fcs_gating_object') {
    print("General note: if this gating tree borrows from another, consider copying the common gate path to this new tree before continuing. Each tree is linear and is independent of all other trees.")
    if(!'gate_trees' %in% names(object)) {
      object[['gate_trees']] <- list()
    }
    if(!tree_name %in% names(object[['gate_trees']])) {
      object[['gate_trees']][[tree_name]] <- list()
    }
    tmpdf <- as.data.frame(object[['data']]); row.names(tmpdf) <- 1:nrow(tmpdf)
    if(parent_name!='none') {
      if(!parent_name %in% names(object[['gate_trees']][[tree_name]])) {
        stop("Could not find 'parent_name' in specified 'tree_name'. If 'parent_name' is not 'none', the specified 'parent_name' must exist in the specified 'tree_name'.")
      } else {
        parent_index <- which(names(object[['gate_trees']][[tree_name]])==parent_name)
        mask_list <- lapply(X = object[['gate_trees']][[tree_name]][1:parent_index], FUN = function(arg) return(arg[['mask']]))
        cell_mask <- do.call(pmin, mask_list)
        # cell_mask <- object[['gate_trees']][[tree_name]][[parent_name]][['mask']]
        df <- tmpdf[cell_mask==1,] # 1 = in parent gate; 0 = not in parent gate
        names(cell_mask) <- row.names(tmpdf)
      }
    } else {
      df <- tmpdf; rm(tmpdf)
    }
  } else if(class(object)=='data.frame') {
    df <- object
  } else if('matrix' %in% class(object)){
    df <- as.data.frame(object)
  } else {
    stop("object should be of class: data.frame, matrix, or fcs_gating_object")
  }
  if(!all(c(a, h) %in% names(df))) {
    stop("df must contain columns '", a, "' and '", h, "'")
  }
  
  x <- df[,a];  y <- df[,h]
  
  # compute the quantile bounds once
  x_lo <- quantile(x, trim_frac)
  x_hi <- quantile(x, 1 - trim_frac)
  y_lo <- quantile(y, trim_frac)
  y_hi <- quantile(y, 1 - trim_frac)
  
  # subset to the central X% in both dims for model fitting
  idx <- which(
    (x >= x_lo & x <= x_hi) &
      (y >= y_lo & y <= y_hi)
  )
  
  # fit the line
  if(length(fit_method)>1) {
    fit_method <- 'rlm'
  }
  fit <- switch(fit_method,
                lm  = lm(y ~ x, data = df[idx,]),
                rlm = rlm(y ~ x, data = df[idx,])
  )
  
  # compute raw residuals
  res <- resid(fit)  # y - ŷ
  rev_res_d <- density(-res)
  # cut the residuals
  if(length(cut_method)>1) {
    cut_method <- 'gmm'
  }
  if(!cut_method %in% c('gmm','flex','both')) {
    stop("cut_method must be one of: 'gmm', 'flex', 'both'")
  }
  cuts <- fcs_get_em_cutpoint(x = -res, general_method = cut_method, curvature_eps = curvature_eps)
  cutpoint <- -cuts$cut

  inside <- res >= cutpoint
  gated  <- df[inside, , drop=FALSE]
  
  fraction_keep <- (sum(inside)/nrow(df))*100
  if(fraction_keep<85) {
    warning(print(paste0('percent of events in singlet gate: ',as.character(round(fraction_keep,1)),'%. Consider adjusting curvature_eps')))
  } else {
    print(paste0('percent of events in singlet gate: ',as.character(round(fraction_keep,1)),'%'))
  }
  if(parent_name!='none') {
    current_gate <- as.numeric(inside); names(current_gate) <- row.names(df)
    pos <- match(names(current_gate), names(cell_mask))
    final_mask <- cell_mask
    final_mask[pos] <- current_gate
    names(final_mask) <- NULL
  } else {
    final_mask <- as.numeric(inside)
  }
  if(class(object)=='fcs_gating_object') {
    new_branch <- list('mask' = final_mask, 
                       'tree' = tree_name, 
                       'parent' = parent_name, 
                       'feature_a' = a, 
                       'feature_h' = h, 
                       'fit_method' = fit_method, 
                       'cut_method' = cut_method, 
                       'gate_fn' = 'fcs_gate_singlets')
    object[['gate_trees']][[tree_name]] <- append(object[['gate_trees']][[tree_name]], list(new_branch))
    names(object[['gate_trees']][[tree_name]])[length(object[['gate_trees']][[tree_name]])] <- gate_name
    return(object)
  } else {
    return(gated)
  }
}

#’ @title Plot Singlet Gate on A vs H
#’
#’ @description
#’   Visualizes area vs height, highlighting singlets in blue and background  
#’   events in black.  
#’
#’ @param object  
#’   An “fcs_gating_object” containing the singlet gate branch.  
#’ @param tree_name  
#’   Character; gating tree name (default “tree1”).  
#’ @param gate_name  
#’   Character; gate branch to plot (default “ssc_singlets”).  
#’ @param alpha  
#’   Numeric; transparency for background points (default 0.25).  
#’ @param psize  
#’   Numeric; point size (default 0.6).  
#’ @param downsample_size  
#’   Integer; maximum points to plot (default 100000).  
#’
#’ @return  
#’   A ggplot2 object.  
#’
#’ @examples
#’ \dontrun{
#’   go <- fcs_gate_singlets(joined, gate_name = "singlets")
#’   p <- fcs_plot_singlets(go, gate_name = "singlets")
#’   print(p)
#’ }
#’
#’ @seealso  
#’   fcs_gate_singlets, ggplot2::ggplot  
#’
#’ @importFrom ggplot2 ggplot aes geom_point theme_bw labs theme element_text  
#’ @export
fcs_plot_singlets <- function(object,
                              tree_name = 'tree1', 
                              gate_name = 'ssc_singlets', 
                              alpha = 0.25,
                              psize = 0.6,
                              downsample_size = 100000) {
  require(ggplot2)
  require(scales)
  
  a <- object[['gate_trees']][[tree_name]][[gate_name]][['feature_a']]
  h <- object[['gate_trees']][[tree_name]][[gate_name]][['feature_h']]
  parent_name <- object[['gate_trees']][[tree_name]][[gate_name]][['parent']]
  cells_all <- object$data[,c(a,h)]
  cells <- cells_all[object[['gate_trees']][[tree_name]][[gate_name]][['mask']]==1,]
  cells_all <- cells_all[object[['gate_trees']][[tree_name]][[parent_name]][['mask']]==1,]
  
  if(!is.na(downsample_size)) {
    ds_scaled <- ceiling((downsample_size/nrow(cells_all))*nrow(cells))
    if(nrow(cells_all)>downsample_size) {
      cells_all <- cells_all[sample(1:nrow(cells_all),downsample_size,replace=F),]
    }
    if(nrow(cells)>downsample_size) {
      cells <- cells[sample(1:nrow(cells),ds_scaled,replace=F),]
    }
  }
  
  capture_cnames <- colnames(cells_all)
  colnames(cells_all) <- c('xcol','ycol'); colnames(cells) <- colnames(cells_all)
  
  pl <- ggplot() + 
    geom_point(data = cells_all, aes(x = xcol, y = ycol), 
               color = 'black', alpha = alpha, size = psize) + 
    geom_point(data = cells, mapping = aes(x = xcol, y = ycol), 
               color = 'blue', alpha = ifelse(alpha*2>1,1,alpha*2), size = psize*1.1) + 
    theme_bw() +
    labs(x = capture_cnames[1], y = capture_cnames[2]) + 
    theme(axis.title = element_text(size = 20, face = 'bold'), 
          axis.text = element_text(size = 14), 
          plot.title = element_text(size = 24, hjust = 0.5, face = 'bold'))
  return(pl)
}

#’ @title EM‐Based Cutpoint Estimation
#’
#’ @description
#’   Fits a two-component Gaussian mixture model (via mclust or mixtools) or  
#’   applies a curvature (“flex”) fallback to estimate a data‐driven cutpoint.  
#’
#’ @param x  
#’   Numeric vector of observations.  
#’ @param em_method  
#’   Character; “mclust” (default) or “mix” for mixtools.  
#’ @param general_method  
#’   Character; “gmm”, “flex”, “mean”, or “both” (default).  
#’ @param post.thresh  
#’   Numeric; posterior probability threshold (default 0.5).  
#’ @param prom_tol  
#’   Numeric; minimum component proportion (default 0.05).  
#’ @param flex_tail  
#’   Numeric; quantile for flex fallback (default 0.75).  
#’ @param bw_adjust  
#’   Numeric; bandwidth adjustment for density (default 2).  
#’ @param n_grid  
#’   Integer; grid length for density (default 512).  
#’ @param curvature_eps  
#’   Numeric; curvature tolerance for flex‐fallback (default 0.01).  
#’ @param return.model  
#’   Logical; if TRUE, include fitted model in output (default TRUE).  
#’
#’ @return  
#’   A list with:  
#’   - cut: numeric cutpoint  
#’   - method: character tag of method used  
#’   - model: the fitted EM object (if return.model = TRUE)  
#’
#’ @examples
#’ \dontrun{
#’   res <- fcs_get_em_cutpoint(x = df$FL1, em_method = "mclust")
#’ }
#’
#’ @seealso  
#’   mclust::Mclust, mixtools::normalmixEM  
#’
#’ @importFrom mclust Mclust emControl  
#’ @importFrom mixtools normalmixEM  
#’ @importFrom stats quantile dnorm density uniroot sd  
#’ @importFrom pracma findpeaks  
#’ @export
fcs_get_em_cutpoint <- function(x, 
                                em_method = c('mclust','mix'),
                                general_method = c('gmm','flex','mean'),
                                post.thresh = 0.5,
                                prom_tol = 0.05,
                                flex_tail = 0.75,
                                bw_adjust = 2,     
                                n_grid = 512,   
                                curvature_eps = 0.01, 
                                return.model = TRUE) {
  em_method <- match.arg(em_method)
  
  # validate general_method
  if (!is.null(general_method) && ! general_method %in% c('gmm','flex','mean')) {
    stop("'general_method' must be one of: 'gmm', 'flex', 'mean', or NULL")
  }
  
  # --- 1) GMM phase (either Mclust or mixtools) ---
  if (is.null(general_method) || general_method %in% c('gmm','mean')) {
    if(is.null(general_method)) {
      general_method <- 'gmm'
    }
    if (em_method=='mclust') {
      # univariate GMM via mclust
      fit <- tryCatch(
        Mclust(x, G = 2, modelNames = 'V',
               control = emControl(itmax = 2000, tol = 1e-8),
               verbose = FALSE),
        error = function(e) NULL
      )
      ok <- !is.null(fit) && all(fit$parameters$pro > prom_tol)
      
      if (ok) {
        pi1 <- fit$parameters$pro[1]
        pi2 <- fit$parameters$pro[2]
        mu1 <- fit$parameters$mean[1]
        mu2 <- fit$parameters$mean[2]
        sd1 <- sqrt(fit$parameters$variance$sigmasq[1])
        sd2 <- sqrt(fit$parameters$variance$sigmasq[2])
      }
      
    } else if (em_method=='mix') {
      # univariate GMM via mixtools
      lower_guess <- quantile(x, 0.05)
      upper_guess <- quantile(x, 0.95)
      fit <- tryCatch(
        normalmixEM(x,
                    k = 2,
                    mu = c(lower_guess, upper_guess),
                    sigma = rep(sd(x),2),
                    lambda = c(0.5,0.5)),
        error = function(e) NULL
      )
      ok <- !is.null(fit) && all(fit$lambda > prom_tol)
      
      if (ok) {
        pi1 <- fit$lambda[1]
        pi2 <- fit$lambda[2]
        mu1 <- fit$mu[1]
        mu2 <- fit$mu[2]
        sd1 <- fit$sigma[1]
        sd2 <- fit$sigma[2]
      }
    }
    
    # if the GMM converged and proportions pass the threshold
    if (ok) {
      λ <- post.thresh/(1-post.thresh)
      h <- function(z) pi1*dnorm(z,mu1,sd1) -
        λ  *pi2*dnorm(z,mu2,sd2)
      lo <- min(mu1, mu2)
      hi <- max(mu1, mu2)
      
      cut_gmm <- tryCatch(
        uniroot(h, lower=lo, upper=hi)$root,
        error = function(e) NA
      )
      
      if (!is.na(cut_gmm)) {
        method_tag <- if (em_method=='mix') 
          sprintf("mix-EM(τ=%.2g)", post.thresh)
        else
          sprintf("Mclust-GMM(τ=%.2g)", post.thresh)
        res <- list(cut = cut_gmm, method = method_tag)
        if (return.model) res$model <- fit
        if(general_method!='mean') {
          return(res)
        }
      }
    }
  }
  ok <- TRUE
  
  # --- 2) Flex‐point fallback ---
  d <- density(x, adjust = bw_adjust, n = n_grid)
  xg <- d$x;   yg <- d$y
  dx <- mean(diff(xg))
  d2 <- c(NA, diff(yg,2)/dx^2, NA)
  absd2 <- abs(d2)
  
  peaks <- findpeaks(-absd2, nups=1, ndowns=1)[,2]
  qt <- quantile(x, flex_tail)
  cand <- peaks[xg[peaks] > qt & absd2[peaks] < curvature_eps]
  
  cut_flex <- if (length(cand)==0) {
    median(x)
  } else {
    xg[cand[which.min(abs(xg[cand] - qt))]]
  }
  
  if(all(ok, general_method=='both')) {
    mean_cut <- mean(c(cut_gmm, cut_flex))
    res <- list(cut = mean_cut, method = 'gmm flex mean', model = fit)
  } else {
    res <- list(cut = cut_flex, method = "flex-fallback")
  }
  return(res)
}

#’ @title Set a Gate Based on EM Cutpoint
#’
#’ @description
#’   Wrapper around fcs_get_em_cutpoint() to apply a univariate gate on  
#’   a single feature, storing positive/negative masks in the gating object.  
#’
#’ @param object  
#’   An “fcs_gating_object” or numeric vector.  
#’ @param feature  
#’   Character; column name of the feature in object$data.  
#’ @param tree_name  
#’   Character; gating tree name (default “tree1”).  
#’ @param parent_name  
#’   Character; parent gate (default “none”).  
#’ @param gate_name  
#’   Character; new gate label (no “+”/“–” allowed).  
#’ @param em_method  
#’   Character; “mclust” or “mix” (default).  
#’ @param general_method  
#’   Character; “gmm”, “flex”, “mean”, or “both” (default).  
#’ @param post.thresh  
#’   Numeric; posterior threshold (default 0.5).  
#’ @param prom_tol  
#’   Numeric; component proportion threshold (default 0.05).  
#’ @param flex_tail  
#’   Numeric; quantile for flex fallback (default 0.75).  
#’ @param bw_adjust  
#’   Numeric; density bandwidth adjust (default 2).  
#’ @param n_grid  
#’   Integer; grid length for density (default 512).  
#’ @param curvature_eps  
#’   Numeric; curvature tolerance (default 0.01).  
#’
#’ @return  
#’   If input is an “fcs_gating_object”, returns it with a new branch under  
#’   `gate_trees[[tree_name]][[gate_name]]` containing positive/negative masks and cutpoint.  
#’   Otherwise returns a list with cutpoint data.  
#’
#’ @seealso  
#’   fcs_get_em_cutpoint  
#’
#’ @export
fcs_set_gate <- function(object,
                         feature, 
                         tree_name = 'tree1', 
                         parent_name = 'cells', 
                         gate_name = 'singlets1', 
                         em_method = c('mclust','mix'),
                         general_method = c('gmm','flex','mean'),
                         post.thresh = 0.5,
                         prom_tol = 0.05,
                         flex_tail = 0.75,
                         bw_adjust = 2,     
                         n_grid = 512,   
                         curvature_eps = 0.01) {
  # fcs_get_em_cutpoint generalized wrapper for class fcs_gating_object
  if(regexpr(pattern = '(\\+|\\-)', text = gate_name)!=-1) {
    stop("'gate_name' should not include '+' or '-'")
  }
  if(class(object)=='fcs_gating_object') {
    print("General note: if this gating tree borrows from another, consider copying the common gate path to this new tree before continuing. Each tree is linear and is independent of all other trees.")
    if(!feature %in% colnames(object[['data']])) {
      stop("Inferred feature must be in column names of 'object[['data']]'")
    }
    if(!'gate_trees' %in% names(object)) {
      object[['gate_trees']] <- list()
    }
    if(!tree_name %in% names(object[['gate_trees']])) {
      object[['gate_trees']][[tree_name]] <- list()
    }
    tmpdf <- object[['data']][,feature]; names(tmpdf) <- 1:length(tmpdf)
    if(parent_name!='none') {
      if(!parent_name %in% names(object[['gate_trees']][[tree_name]])) {
        stop("Could not find 'parent_name' in specified 'tree_name'. If 'parent_name' is not 'none', the specified 'parent_name' must exist in the specified 'tree_name'.")
      } else {
        parent_index <- which(names(object[['gate_trees']][[tree_name]])==parent_name)
        mask_list <- lapply(X = object[['gate_trees']][[tree_name]][1:parent_index], FUN = function(arg) return(arg[['mask']]))
        cell_mask <- do.call(pmin, mask_list); names(cell_mask) <- names(tmpdf)
        df <- tmpdf[cell_mask==1]
      }
    } else {
      df <- tmpdf; rm(tmpdf)
    }
  } else if(class=='numeric'){
    df <- object
  } else {
    stop("object should be of class: numeric (length N, where N>>1) or fcs_gating_object")
  }
  if(length(em_method)>1) {
    em_method <- 'mix'
  }
  if(length(general_method)>1) {
    general_method <- 'gmm'
  }
  # validate general_method
  if (!is.null(general_method) && ! general_method %in% c('gmm','flex','mean')) {
    stop("'general_method' must be one of: 'gmm', 'flex', 'mean', or NULL")
  }
  # GMM phase (either mclust::Mclust or mixtools::normalmixEM)
  if (is.null(general_method) || general_method %in% c('gmm','mean')) {
    if (em_method=='mclust') {
      # univariate GMM via mclust
      fit <- tryCatch(
        mclust::Mclust(df, G = 2, modelNames = 'V',
               control = emControl(itmax = 2000, tol = 1e-8),
               verbose = FALSE),
        error = function(e) NULL
      )
      ok <- !is.null(fit) && all(fit$parameters$pro > prom_tol)
      if (ok) {
        pi1 <- fit$parameters$pro[1]
        pi2 <- fit$parameters$pro[2]
        mu1 <- fit$parameters$mean[1]
        mu2 <- fit$parameters$mean[2]
        sd1 <- sqrt(fit$parameters$variance$sigmasq[1])
        sd2 <- sqrt(fit$parameters$variance$sigmasq[2])
      }
      
    } else if (em_method=='mix') {
      # univariate GMM via mixtools
      # lower_guess <- quantile(df, 0.05)
      # upper_guess <- quantile(df, 0.95)
      fit <- tryCatch(
        mixtools::normalmixEM(df,
                    k = 2,
                    # mu = c(lower_guess, upper_guess),
                    sigma = rep(sd(df),2),
                    lambda = c(0.5,0.5)),
        error = function(e) NULL
      )
      ok <- !is.null(fit) && all(fit$lambda > prom_tol)
      if (ok) {
        pi1 <- fit$lambda[1]
        pi2 <- fit$lambda[2]
        mu1 <- fit$mu[1]
        mu2 <- fit$mu[2]
        sd1 <- fit$sigma[1]
        sd2 <- fit$sigma[2]
      }
    }
    # if the GMM converged and proportions pass the threshold
    if (ok) {
      λ <- post.thresh/(1-post.thresh)
      h <- function(z) pi1*dnorm(z,mu1,sd1) -
        λ  *pi2*dnorm(z,mu2,sd2)
      lo <- min(mu1, mu2)
      hi <- max(mu1, mu2)
      
      cut_gmm <- tryCatch(
        uniroot(h, lower=lo, upper=hi)$root,
        error = function(e) NA
      )
      
      if (!is.na(cut_gmm)) {
        method_tag <- if (em_method=='mix') {
          sprintf("mix-EM(τ=%.2g)", post.thresh)
        } else {
          sprintf("Mclust-GMM(τ=%.2g)", post.thresh)
        }
        res <- list(cut = cut_gmm, method = method_tag, model = fit)
        if(general_method!='mean') {
          if(class(object)=='fcs_gating_object') {
            cutpt <- cut_gmm
          } else {
            return(res)
          }
        }
      }
    }
  } else {
    ok <- FALSE
  }
  if(any(!ok, general_method=='mean')) {
    # flex‐point fallback
    d <- density(df, adjust = bw_adjust, n = n_grid)
    xg <- d$x; yg <- d$y
    dx <- mean(diff(xg))
    d2 <- c(NA, diff(yg,2)/dx^2, NA)
    absd2 <- abs(d2)
    
    peaks <- pracma::findpeaks(-absd2, nups=1, ndowns=1)[,2]
    qt <- quantile(df, flex_tail)
    cand <- peaks[xg[peaks] > qt & absd2[peaks] < curvature_eps]
    
    cut_flex <- if (length(cand)==0) {
      median(df)
    } else {
      xg[cand[which.min(abs(xg[cand] - qt))]]
    }
    
    if(all(ok, general_method=='mean')) {
      cutpt <- mean(c(cut_gmm, cut_flex))
    } else {
      cutpt <- cut_flex
    }
  }
  if(class(object)=='fcs_gating_object') {
    inside_positive <- df >= cutpt
    inside_negative <- df < cutpt
    if(parent_name!='none') {
      # positive
      current_gate_pos <- as.numeric(inside_positive); names(current_gate_pos) <- names(df)
      pos <- match(names(current_gate_pos), names(cell_mask))
      positive_mask <- cell_mask
      positive_mask[pos] <- current_gate_pos
      names(positive_mask) <- NULL
      # negative
      current_gate_negative <- as.numeric(inside_negative); names(current_gate_negative) <- names(df)
      pos <- match(names(current_gate_negative), names(cell_mask))
      negative_mask <- cell_mask
      negative_mask[pos] <- current_gate_negative
      names(negative_mask) <- NULL
    } else {
      positive_mask <- as.numeric(inside_positive)
      negative_mask <- as.numeric(inside_negative)
    }
    new_branch <- list('positive_mask' = positive_mask, 
                       'negative_mask' = negative_mask,
                       'threshold' = cutpt, 
                       'tree' = tree_name, 
                       'parent' = parent_name, 
                       'feature' = feature,  
                       'em_method' = em_method,
                       'general_method' = general_method,
                       'post.thresh' = post.thresh,
                       'prom_tol' = prom_tol,
                       'flex_tail' = flex_tail,
                       'bw_adjust' = bw_adjust,     
                       'n_grid' = n_grid, 
                       'curvature_eps' = curvature_eps,
                       'gate_fn' = 'fcs_set_gate')
    if(ok) {
      new_branch[['fit']] <- fit
    } else {
      new_branch[['fit']] <- 'none'
    }
    object[['gate_trees']][[tree_name]] <- append(object[['gate_trees']][[tree_name]], list(new_branch))
    names(object[['gate_trees']][[tree_name]])[length(object[['gate_trees']][[tree_name]])] <- gate_name
    return(object)
  } else {
    if(all(ok, general_method=='both')) {
      res <- list(cut = cutpt, method = 'gmm flex mean', model = fit)
    } else {
      if(ok) {
        res <- list(cut = cutpt, method = "flex-fallback", model = fit)
      } else {
        res <- list(cut = cutpt, method = "flex-fallback", model = 'none')
      }
    }
  }
}

#’ @title Define a Phenotypic Node by Boolean Gates
#’
#’ @description
#’   Combines existing “+”/“–” gates in a tree into a composite phenotype node,  
#’   storing the logical AND of positive/negative masks.  
#’
#’ @param object  
#’   An “fcs_gating_object” with completed gate_trees.  
#’ @param phenotype  
#’   Character; expression like “A+B–C+…” where A,B,C are gate names.  
#’ @param tree_name  
#’   Character; gating tree name (default “tree1”).  
#’ @param node_name  
#’   Character; new node label (no “+”/“–” allowed).  
#’
#’ @return  
#’   The updated “fcs_gating_object” with `gate_trees[[tree_name]][[node_name]]`  
#’   containing the composite mask and phenotype metadata.  
#’
#’ @seealso  
#’   fcs_set_gate  
#’
#’ @export
fcs_set_node <- function(object, 
                         phenotype, 
                         tree_name = 'tree1', 
                         node_name = 'phenotype1') {
  # phenotype may be given as A+B-C+D+E- etc where perceived features (A, B, C, D, E) are names of branches within 'tree_name'
  if(regexpr(pattern = '(\\+|\\-)', text = node_name)!=-1) {
    stop("'node_name' should not include '+' or '-'... 'phenotype' is stored within the branch data.")
  }
  get_direction <- strsplit(x = phenotype, split = '[A-Za-z0-9\\.]+')[[1]]
  get_direction <- get_direction[!get_direction=='']
  get_gates <- strsplit(x = phenotype, split = '(\\+|\\-)')[[1]]
  gate_dir_df <- data.frame(name = get_gates, direction = get_direction)
  
  if(mean(get_gates %in% names(object[['gate_trees']][[tree_name]]))!=1) {
    print(paste0("Features not found: ",paste0(get_gates[!get_gates %in% names(object[['gate_trees']][[tree_name]])], collapse=",")))
    stop("All features given in 'phenotype' must be present in names(object[['gate_trees']][[tree_name]]).")
  }
  mask_list <- list()
  for(i in 1:nrow(gate_dir_df)) {
    if(gate_dir_df$direction[i]=='+') {
      mask_list[[i]] <- object$gate_trees[[tree_name]][[gate_dir_df$name[i]]]$positive_mask
    } else if(gate_dir_df$direction[i]=='-') {
      mask_list[[i]] <- object$gate_trees[[tree_name]][[gate_dir_df$name[i]]]$negative_mask
    }
  }
  data_mask <- do.call(pmin, mask_list) # 1 = in node; 0 = not in node
  
  new_branch <- list('mask' = data_mask, 
                     'tree' = tree_name,  
                     'phenotype' = phenotype, 
                     'feature' = paste0(get_gates,collapse=','), 
                     'direction' = paste0(get_direction,collapse=','), 
                     'gate_fn' = 'fcs_set_node')
  object[['gate_trees']][[tree_name]] <- append(object[['gate_trees']][[tree_name]], list(new_branch))
  names(object[['gate_trees']][[tree_name]])[length(object[['gate_trees']][[tree_name]])] <- node_name
  return(object)
}

#’ @title Plot GMM Fit and Cutpoint
#’
#’ @description
#’   Displays the kernel density of x, overlays the two GMM components and mixture  
#’   density, and marks the chosen cutpoint.  
#’
#’ @param x  
#’   Numeric vector used to fit the GMM.  
#’ @param fit  
#’   Fitted Mclust or normalmixEM model.  
#’ @param cut  
#’   Numeric; cutpoint to highlight (optional).  
#’ @param linewidth  
#’   Numeric; line width for plotting (default 3).  
#’
#’ @return  
#’   Invisibly `NULL`; produces a base‐graphics plot.  
#’
#’ @examples
#’ \dontrun{
#’   fit <- mclust::Mclust(x, G=2)
#’   fcs_plot_gmm_fit(x, fit, cut=1.2)
#’ }
#’
#’ @seealso  
#’   fcs_get_em_cutpoint, stats::density  
#’
#’ @importFrom stats density dnorm  
#’ @export
fcs_plot_gmm_fit <- function(x, fit, cut = NULL, linewidth = 3) {
  d <- density(x, n = 512)
  plot(d, main = "Data + Fitted GMM", lwd = linewidth, col = "black")
  
  if (!is.null(fit$parameters)) {
    # an Mclust object
    p  <- fit$parameters$pro
    mu <- fit$parameters$mean
    sd <- sqrt(fit$parameters$variance$sigmasq)
  } else if (! is.null(fit$lambda)) {
    # a mixtools::normalmixEM object
    p  <- fit$lambda
    mu <- fit$mu
    # note: fit$sigma is already the standard deviations
    sd <- fit$sigma
  } else {
    stop("plot_gmm: unrecognized model object. ",
         "Pass an Mclust or normalmixEM fit.")
  }
  
  xs <- seq(min(x), max(x), length = 500)
  y1 <- p[1] * dnorm(xs, mu[1], sd[1])
  y2 <- p[2] * dnorm(xs, mu[2], sd[2])
  lines(xs, y1,      col = "firebrick", lwd = linewidth, lty = 2)
  lines(xs, y2,      col = "steelblue", lwd = linewidth, lty = 2)
  lines(xs, y1 + y2, col = "darkgreen", lwd = linewidth, lty = 3)
  
  if (!is.null(cut)) {
    abline(v = cut, col = "purple", lwd = linewidth, lty = 4)
    legend("topright",
           legend = c("KDE", "Comp1", "Comp2", "Mixture", "Cut"),
           col    = c("black","firebrick","steelblue","darkgreen","purple"),
           lty    = c(1,2,2,3,4),
           lwd    = rep(linewidth,5),
           bty    = "n")
  } else {
    legend("topright",
           legend = c("KDE", "Comp1", "Comp2", "Mixture"),
           col    = c("black","firebrick","steelblue","darkgreen"),
           lty    = c(1,2,2,3),
           lwd    = rep(linewidth,4),
           bty    = "n")
  }
}

#’ @title Plot Quadrant Percentages on Two‐Dimensional Data
#’
#’ @description
#’   Splits data.frame `df` into four quadrants by `xcut`/`ycut`, computes  
#’   percentages in UL/UR/LR/LL, and overlays density contours and labels.  
#’
#’ @param df  
#’   data.frame with two numeric columns to be split.  
#’ @param xcut  
#’   Numeric; vertical cut‐point on first column.  
#’ @param ycut  
#’   Numeric; horizontal cut‐point on second column.  
#’ @param alpha  
#’   Numeric; transparency for points (default 0.25).  
#’ @param psize  
#’   Numeric; point size (default 0.6).  
#’ @param n  
#’   Integer; grid size for density (default 100).  
#’ @param bins  
#’   Integer; contour bin count (default 12).  
#’
#’ @return  
#’   A ggplot2 object with points, contours, quadrant lines, and percentage labels.  
#’
#’ @examples
#’ \dontrun{
#’   df <- data.frame(x=rnorm(1e4), y=rnorm(1e4))
#’   p <- fcs_plot_quadrants(df, xcut=0, ycut=0)
#’   print(p)
#’ }
#’
#’ @seealso  
#’   ggplot2::stat_density_2d  
#’
#’ @importFrom ggplot2 ggplot aes geom_point stat_density_2d geom_vline geom_hline geom_text theme_bw labs theme element_text  
#’ @export
fcs_plot_quadrants <- function(df,
                               xcut, 
                               ycut,
                               alpha = 0.25,
                               psize = 0.6,
                               n = 100,
                               bins = 12) {
  require(ggplot2)
  require(scales)
  
  downsample_size <- 100000
  if(nrow(df)>downsample_size) {
    set.seed(123); df <- df[sample(1:nrow(df),downsample_size,replace=F),]
  }
  x <- df[,1];  y <- df[,2]
  tot <- length(x)
  pct_ul <- sum(x < xcut & y >= ycut)/tot * 100
  pct_ur <- sum(x >= xcut & y >= ycut)/tot * 100
  pct_lr <- sum(x >= xcut & y <  ycut)/tot * 100
  pct_ll <- sum(x < xcut & y <  ycut)/tot * 100
  labs   <- sprintf("%.1f%%", c(pct_ul, pct_ur, pct_lr, pct_ll))
  
  pos_df <- data.frame(x = c(min(x), max(x), max(x), min(x)),
                       y = c(max(y), max(y), min(y), min(y)),
                       label = paste0(paste0('Q',1:4,'\n'),labs),
                       hjust = c(0, 1, 1, 0),
                       vjust = c(1, 1, 0, 0))
  capture_cnames <- colnames(df)
  colnames(df) <- c('xcol','ycol')
  
  h <- c(diff(range(df[,1]))/30,
         diff(range(df[,2]))/30)
  
  pl <- ggplot(df, aes(x = xcol, y = ycol)) + 
    geom_point(color = "black", alpha = alpha, size = psize) + 
    stat_density_2d(geom = "contour",
                    color = "cyan",
                    h = h,
                    n = n,
                    bins = bins,
                    size = 0.75) + 
    geom_vline(xintercept = xcut,
               linetype = 6,
               color = "purple",
               size = 1.2) +
    geom_hline(yintercept = ycut,
               linetype = 6,
               color = "purple",
               size = 1.2) + 
    geom_text(data = pos_df,
              aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE,
              size = 6,
              fontface = 'bold', 
              color  = "black") + 
    theme_bw() +
    labs(x = capture_cnames[1], y = capture_cnames[2]) + 
    theme(axis.title = element_text(size = 20, face = 'bold'), 
          axis.text = element_text(size = 14), 
          plot.title = element_text(size = 24, hjust = 0.5, face = 'bold'))
  return(pl)
}
