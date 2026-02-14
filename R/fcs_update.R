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
#' @details
#'   - If both `collection_instrument` and `object_history` are present in
#'     `fcs_join_obj`, a warning `"object is already up to date"` is issued
#'     and the original object is returned.
#'   - Otherwise, `collection_instrument` is auto-detected based on the data:
#'     if min(data) >= 0, assumes 'cytof' (mass cytometry); otherwise assumes
#'     'flow' (fluorescence cytometry). `object_history` is initialized with
#'     `"updated: <timestamp>"`.
#'
#' @return
#'   The input `fcs_join_obj`, with the fields:
#'   - `collection_instrument`: auto-detected as 'cytof' or 'flow' based on data.
#'   - `object_history`: a character string recording the update timestamp.
#'
#' @examples
#' \dontrun{
#'   files <- list.files("data/fcs", "\\.fcs$", full.names = TRUE)
#'   joined <- FCSimple::fcs_join(files)
#'
#'   # Update adds metadata (instrument_type is auto-detected)
#'   updated <- FCSimple::fcs_update(joined)
#'
#'   # Calling again warns and returns unchanged
#'   unchanged <- FCSimple::fcs_update(updated)
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, FCSimple::fcs_gating_object
#'
#' @export
fcs_update <- function(fcs_join_obj)
{
  if(system(command = 'python --version')==0) {
    pyv <- system(command = 'python --version', intern = TRUE)
    # assumes pip is installed if python is installed; pip3 fallback if pip fails
    pipl <- tryCatch(system("pip list", intern = TRUE), error = function(e) character())
    if (length(pipl) == 0) {
      # try pip3 if pip unsuccessful
      pipl <- tryCatch(system("pip3 list", intern = TRUE), error = function(e) character())
      pipl <- gsub(pattern = '( )+', replacement = 'UniqueSeparator', x = pipl)[-c(1,2)]
      pipl_spl <- strsplit(x = pipl, split = 'UniqueSeparator')
      pip_df <- data.frame(Package = sapply(pipl_spl,function(x) return(x[1])), 
                           Version = sapply(pipl_spl,function(x) return(x[2])))
    } else {
      pipl <- pipl[-c(1,2)]
      split_lines <- strsplit(pipl, "\\s+")
      pip_df <- do.call(rbind, lapply(split_lines, function(x) x[1:2]))
      pip_df <- as.data.frame(pip_df, stringsAsFactors = FALSE)
      colnames(pip_df) <- c("Package", "Version")
    }
  } else {
    pyv <- 'none'
    pip_df <- 'none'
  }
  rv <- getRversion()
  fcs_join_obj[['versions']] <- list(R = rv, 
                                     Python = pyv, 
                                     Rsession = sessionInfo(), 
                                     pip_list = pip_df)
  if(all('object_history' %in% names(fcs_join_obj),
         "collection_instrument" %in% names(fcs_join_obj),
         'run_date' %in%  names(fcs_join_obj),
         'metadata' %in% names(fcs_join_obj))) { # check if up to date first
    warning("object is already up to date")
    return(fcs_join_obj)
  } else {
    if(mean(c("data", "source", "run_date") %in% names(fcs_join_obj))!=1) {
      stop('fcs_join_obj must contain: "data", "source", and "run_date"')
    }
    if(!"metadata" %in% names(fcs_join_obj)) {
      src_data <- fcs_join_obj[['source']]
      base_metadata <- data.frame(patient_ID = gsub(pattern = '_[0-9]+\\-[A-Za-z]+\\-[0-9]+|(\\.csv$|\\.fcs$)',
                                                    replacement = '',
                                                    x = src_data),
                                  run_date = fcs_join_obj[['run_date']])
      fcs_join_obj[['metadata']] <- base_metadata[!duplicated(base_metadata$patient_ID),]
    }
    
    # Auto-detect instrument type based on data
    detected_instrument <- if(min(fcs_join_obj$data, na.rm = TRUE) >= 0) "cytof" else "flow"
    
    fcs_join_obj[['collection_instrument']] <- detected_instrument
    fcs_join_obj[['object_history']] <- paste0("updated: ",Sys.time())
    return(fcs_join_obj)
  }
}
