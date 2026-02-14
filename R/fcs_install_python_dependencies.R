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
#' @param python Character; path or command name for the Python interpreter.
#' @return Logical; TRUE if all required packages are present (or installed),
#'   FALSE otherwise. Invisibly returns result and prints summary messages.
#' @export
fcs_install_python_dependencies <- function(install = FALSE, precompile = FALSE, python = "python") {
  script <- system.file("python", "fcs_install_python_deps.py", package = "FCSimple")
  if (script == "") stop("Installer script not found in package inst/python")
  args <- character()
  if (install) args <- c(args, "--install")
  if (precompile) args <- c(args, "--precompile")
  cmd <- c(shQuote(script), args)
  res <- tryCatch({
    out <- system2(command = python, args = cmd, stdout = TRUE, stderr = TRUE)
    cat(paste(out, collapse = "\n"), "\n")
    TRUE
  }, error = function(e) {
    message("Error running installer: ", e$message)
    FALSE
  })
  invisible(res)
}
