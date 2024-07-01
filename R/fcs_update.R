fcs_update <- function(fcs_join_obj)
{
  if('object_history' %in% names(fcs_join_obj)) { # check if up to date first
    stop("object is already up to date")
  } else {
    fcs_join_obj[['object_history']] <- paste0("updated: ",Sys.time())
    return(fcs_join_obj)
  }
}
