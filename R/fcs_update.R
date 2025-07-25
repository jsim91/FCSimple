#' @title Update FCSimple Object with Instrument Metadata
#'
#' @description
#'   Adds or updates the collection instrument type and history in an FCSimple  
#'   analysis object. If the object already contains both `collection_instrument`  
#'   and `object_history`, the function issues a warning and returns the object  
#'   unchanged.
#'
#' @param fcs_join_obj
#'   An FCSimple analysis object (list) returned by `FCSimple::fcs_join()`.  
#'   May already include `collection_instrument` and `object_history`.
#'
#' @param instrument_type
#'   Character of length one; the type of cytometry data. Must be either  
#'   `"cytof"` (mass cytometry) or `"flow"` (fluorescence cytometry).  
#'   Default is the first element of `c("cytof","flow")`.
#'
#' @details
#'   - If both `collection_instrument` and `object_history` are present in  
#'     `fcs_join_obj`, a warning `"object is already up to date"` is issued  
#'     and the original object is returned.  
#'   - Otherwise, `collection_instrument` is set to `instrument_type`, and  
#'     `object_history` is initialized with `"updated: <timestamp>"`.
#'
#' @return
#'   The input `fcs_join_obj`, with the fields:
#'   - `collection_instrument`: set to the specified `instrument_type`.
#'   - `object_history`: a character string recording the update timestamp.
#'
#' @examples
#' \dontrun{
#'   files <- list.files("data/fcs", "\\.fcs$", full.names = TRUE)
#'   joined <- FCSimple::fcs_join(files)
#'
#'   # First update adds metadata
#'   updated <- FCSimple::fcs_update(joined, instrument_type = "flow")
#'
#'   # Calling again warns and returns unchanged
#'   unchanged <- FCSimple::fcs_update(updated, instrument_type = "cytof")
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, FCSimple::fcs_gating_object
#'
#' @export
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
