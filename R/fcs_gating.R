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

#’ @title    EM Cut‐point Estimation
#’ @description
#’ Fit a two‐component GMM via either mclust or mixtools and find the posterior cutpoint.
#’
#’ @param x              Numeric vector of observations.
#’ @param em_method      Which EM backend: 'mclust' or 'mix'.
#’ @param general_method GMM‐based ('gmm') or flex‐point fallback ('flex').
#’ @param post.thresh    Posterior threshold for cut‐point (default 0.5).
#’ @param prom_tol       Minimum component proportion to accept (default 0.05).
#’ @param flex_tail      Quantile for flex‐point fallback (default 0.75).
#’ @param bw_adjust      Bandwidth adjustment for density() (default 2).
#’ @param n_grid         Grid length for density() (default 512).
#’ @param curvature_eps  Curvature tolerance for flex‐point (default 0.01).
#’ @param return.model   If TRUE, return the fitted EM model in the output.
#’
#’ @return A list with components  
#’   * cut: the numeric cut‐point  
#’   * method: character tag  
#’   * model: (optional) the EM fit object  
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
