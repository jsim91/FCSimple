#’ @title Hyperlog Transform Raw Cytometry Data
#’
#’ @description
#’   Applies a hyperlog transformation to the raw cytometry measurements in an
#’   FCSimple joined object. Uses flowCore::hyperlogtGml2 internally to convert
#’   linear raw data into a hyperlog scale that accommodates both low and high
#’   intensity signals.
#’
#’ @param fcs_join_obj
#’   A list returned by FCSimple::fcs_join(), containing at minimum:
#’   - raw_data: numeric matrix or data.frame of events × channels  
#’   - source: character vector indicating sample origin for each event
#’
#’ @param hyperlog_transform_T
#’   Numeric; top‐of‐scale parameter T controlling the maximum display value
#’   (default 100000).
#’
#’ @param hyperlog_transform_M
#’   Numeric; number of decades to transform (default 5).
#’
#’ @param hyperlog_transform_W
#’   Numeric; width of the linear region around zero (default 0.001).
#’
#’ @param hyperlog_transform_A
#’   Numeric; additional constant to shift the transform (default 2).
#’
#’ @details
#’   The function splits `raw_data` by `source`, applies the hyperlogtGml2
#’   transformation to each channel within each sample, and then recombines
#’   the transformed data in the original order. The hyperlog parameters
#’   (T, M, W, A) control the shape of the scale, with larger T/M extending
#’   the dynamic range and W/A tuning the linear-to-log transition.
#’
#’ @return
#’   The input `fcs_join_obj`, with its `data` element replaced by the
#’   hyperlog‐transformed matrix (dimensions identical to `raw_data`).
#’
#’ @examples
#’ \dontrun{
#’   joined <- FCSimple::fcs_join(list(ff1, ff2))
#’   hyper_obj <- FCSimple::fcs_as.hyperlog(
#’     joined,
#’     hyperlog_transform_T = 1e5,
#’     hyperlog_transform_M = 5,
#’     hyperlog_transform_W = 0.001,
#’     hyperlog_transform_A = 2
#’   )
#’ }
#’
#’ @seealso
#’   flowCore::hyperlogtGml2, FCSimple::fcs_join
#’
#’ @importFrom flowCore hyperlogtGml2
#’ @export
fcs_as.hyperlog <- function(fcs_join_obj,
                            hyperlog_transform_T = 100000,
                            hyperlog_transform_M = 5,
                            hyperlog_transform_W = 0.001,
                            hyperlog_transform_A = 2) {
  # Split raw_data by sample
  myfs <- split(
    x = as.data.frame(fcs_join_obj$raw_data),
    f = fcs_join_obj$source
  )
  
  # Internal helper: apply hyperlog transform to each sample
  transf_hyperlog <- function(fset, hyper_t, hyper_m, hyper_w, hyper_a) {
    for (i in seq_along(fset)) {
      exprs_data <- fset[[i]]
      for (j in seq_len(ncol(exprs_data))) {
        # build the hyperlog transform function for this channel
        transform_fun <- flowCore::hyperlogtGml2(
          parameters     = colnames(exprs_data)[j],
          'T'            = hyper_t,
          M              = hyper_m,
          W              = hyper_w,
          A              = hyper_a,
          transformationId = "hyper1"
        )
        # apply transform
        exprs_data[, j] <- eval(transform_fun)(exprs_data)
      }
      # accumulate
      if (i == 1) {
        transf_data <- exprs_data
      } else {
        transf_data <- rbind(transf_data, exprs_data)
      }
    }
    transf_data
  }
  
  # perform the transformation
  tmp_data <- transf_hyperlog(
    fset    = myfs,
    hyper_t = hyperlog_transform_T,
    hyper_m = hyperlog_transform_M,
    hyper_w = hyperlog_transform_W,
    hyper_a = hyperlog_transform_A
  )
  
  # replace data slot and return
  fcs_join_obj$data <- tmp_data
  return(fcs_join_obj)
}
