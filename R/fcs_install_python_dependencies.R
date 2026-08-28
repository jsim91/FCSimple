#' Install or check Python dependencies for FCSimple
#'
#' This function calls the bundled Python helper `fcs_install_python_deps.py`
#' which checks for required Python packages and can optionally install them
#' using the same Python interpreter that is invoked. The helper avoids
#' importing heavy native libraries to prevent JIT/compile side effects.
#'
#' @param install Logical; if TRUE, attempt to install missing packages via pip.
#' @param precompile Logical; if TRUE, run a one-time single-threaded precompile
#'   of numba/umap to warm the JIT cache (optional).
#' @param build_optsne Logical; if TRUE, build and install the bundled opt-SNE
#'   (MulticoreTSNE) package from source. Requires cmake and a C++ compiler.
#' @param python Character; path or command name for the Python interpreter.
#' @return Logical; TRUE if all required packages are present (or installed),
#'   FALSE otherwise. Invisibly returns result and prints summary messages.
#' @export
fcs_install_python_dependencies <- function(install = FALSE, precompile = FALSE,
                                            build_optsne = FALSE, python = "python") {
  script <- system.file("python", "fcs_install_python_deps.py", package = "FCSimple")
  if (script == "") stop("Installer script not found in package inst/python")
  args <- character()
  if (install) args <- c(args, "--install")
  if (precompile) args <- c(args, "--precompile")
  if (build_optsne) args <- c(args, "--build-optsne")
  cmd <- c(shQuote(script), args)
  res <- tryCatch({
    # system2 emits a warning when the command exits non-zero; the Python
    # helper exits 1 to report missing packages (when install=FALSE), which is
    # an expected outcome rather than an R error. Suppress that warning and
    # translate the exit status into the documented return value below.
    out <- suppressWarnings(system2(command = python, args = cmd, stdout = TRUE, stderr = TRUE))
    cat(paste(out, collapse = "\n"), "\n")
    status <- attr(out, "status")
    ok <- is.null(status) || isTRUE(status == 0)
    if (!ok) {
      message("FCSimple Python setup reported an error (see output above).")
    }
    ok
  }, error = function(e) {
    message("Error running installer: ", e$message)
    FALSE
  })
  invisible(res)
}
