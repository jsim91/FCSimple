# Internal helpers for locating writable temporary directories.
#
# Python-backed steps (clustering, dimension reduction, hyperlog transforms,
# and the per-channel transform app) exchange scratch files with Python. These
# must live in a writable location because the installed package directory is
# typically read-only (and may be shared), so each call gets a unique directory
# under tempdir() that is removed afterwards.

#' Create a unique temporary directory for interop scratch files.
#'
#' @return Character path to a directory that is guaranteed to exist.
#' @keywords internal
#' @noRd
.fcs_temp_dir <- function() {
  d <- tempfile("FCSimple_")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Return the feature (column) names of a feature matrix.
#'
#' @param m A matrix or data frame of events x features.
#' @return Character vector of feature names; synthesises names when the
#'   matrix is unnamed.
#' @keywords internal
#' @noRd
.fcs_features <- function(m) {
  cn <- colnames(m)
  if (is.null(cn)) cn <- paste0("feature_", seq_len(ncol(m)))
  cn
}
