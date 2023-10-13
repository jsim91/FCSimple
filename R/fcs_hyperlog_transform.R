fcs_as.hyperlog <- function(fcs_join_obj, hyperlog_transform_T = 100000,
                            hyperlog_transform_M = 5,
                            hyperlog_transform_W = 0.001,
                            hyperlog_transform_A = 2) {

  fset <- split(x = fcs_join_obj$raw_data, f = fcs_join_obj$source)

  transf_hyperlog <- function(fset, hyper_t, hyper_m, hyper_w, hyper_a) {
    for(i in 1:length(fset)) {
      exprs_data <- fset[[i]]
      for(j in 1:ncol(exprs_data)) {
        transform_fun <- flowCore::hyperlogtGml2(parameters = colnames(exprs_data)[j], 'T' = hyper_t, M = hyper_m, W = hyper_w, A = hyper_a, transformationId = "hyper1")
        exprs_data[,j] <- eval(transform_fun)(exprs_data)
      }
      if(i==1) {
        transf_data <- exprs_data
      } else {
        transf_data <- rbind(transf_data, exprs_data)
      }
    }
    return(transf_data)
  }

  tmp_data <- transf_hyperlog(fset = fs, hyper_t = hyperlog_transform_T, hyper_m = hyperlog_transform_M,
                              hyper_w = hyperlog_transform_W, hyper_a = hyperlog_transform_A)

  fcs_join_obj$data <- tmp_data
  return(fcs_join_obj)
}
