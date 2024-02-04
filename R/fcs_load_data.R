fcs_load_data <- function() {
  package_dir <- gsub("\\/","",system.file(package = "FCSimple"))
  data_loc <- paste0(package_dir,"/data")
  fcs_fil <- list.files(path = data_loc, pattern = "fcs$", full.names = TRUE, include.dirs = TRUE)
  return(fcs_fil)
}
