fcs_update <- function(fcs_join_obj, instrument_type = c("cytof","flow"))
{
  if(all('object_history' %in% names(fcs_join_obj), "collection_instrument" %in% names(fcs_join_obj))) { # check if up to date first
    warning("object is already up to date")
    return(fcs_join_obj)
  } else {
    if(any(length(instrument_type)!=1, !instrument_type %in% c("cytof","flow"))) {
      stop("'collection_instrument should be one of: 'cytof' for mass cytometry or 'flow' for fluorescence cytometry")
    } else {
      fcs_join_obj[['collection_instrument']] <- instrument_type
      fcs_join_obj[['object_history']] <- paste0("updated: ",Sys.time())
      return(fcs_join_obj)
    }
  }
}
