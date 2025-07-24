#' @title Load Example FCS Files Bundled with FCSimple
#'
#' @description
#'   Retrieves the file paths of all .fcs files included in the FCSimple
#'   package’s internal `data/` directory. Useful for quickly accessing
#'   example or test files distributed with the package.
#'
#' @return
#'   A character vector of full file paths to each .fcs file found in
#'   `system.file("data", package = "FCSimple")`.
#'
#' @examples
#' \dontrun{
#' # List all example FCS files shipped with FCSimple
#' files <- FCSimple::fcs_load_data()
#'
#' # Use them to join and preprocess
#' joined <- FCSimple::fcs_join(files)
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, utils::system.file
#'
#' @export
fcs_load_data <- function() {
  package_dir <- gsub("\\/$","",system.file(package = "FCSimple"))
  data_loc <- paste0(package_dir,"/data")
  fcs_fil <- list.files(path = data_loc, pattern = "fcs$", full.names = TRUE, include.dirs = TRUE)
  return(fcs_fil)
}
