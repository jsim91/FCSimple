fcs_update <- function(fcs_join_obj, collection_instrument = instrument_type)
{
  if(all('object_history' %in% names(fcs_join_obj), "collection_instrument" %in% names(fcs_join_obj))) { # check if up to date first
    stop("object is already up to date")
    return(fcs_join_obj)
  } else {
    if(any(length(collection_instrument)!=1, !collection_instrument %in% c("cytof","flow"))) {
      stop("'collection_instrument should be one of: 'cytof' or 'flow'")
    } else {
      fcs_join_obj[['collection_instrument']] <- collection_instrument
      fcs_join_obj[['object_history']] <- paste0("updated: ",Sys.time())
      return(fcs_join_obj)
    }
  }
}
