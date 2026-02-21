#' @title Update an FCSimple Object with Metadata, Instrument Type, and Version Info
#'
#' @description
#'   Finalises an FCSimple analysis object produced by \code{FCSimple::fcs_join()}
#'   by performing four tasks:
#'   \enumerate{
#'     \item Captures the current R version, R session info, and (if available)
#'           Python version and installed pip packages, storing them in
#'           \code{obj$versions}.
#'     \item Validates that the object contains the required fields
#'           \code{"data"}, \code{"source"}, and \code{"run_date"}.
#'     \item Validates or constructs per-sample metadata. If the object already
#'           contains a \code{metadata} \code{data.frame} with a \code{patient_ID}
#'           column that fully and exactly covers all unique values in
#'           \code{source} (no missing, no extra, no duplicates), the existing
#'           metadata is preserved. If \code{run_date} is absent from the
#'           existing metadata it is joined in from \code{source}/\code{run_date}.
#'           If any validation check fails, metadata is rebuilt from scratch as a
#'           deduplicated \code{data.frame} with columns \code{patient_ID} and
#'           \code{run_date}.
#'     \item Auto-detects the collection instrument type: if
#'           \code{min(obj$data) >= 0} the instrument is assumed to be
#'           \code{"cytof"} (mass cytometry); otherwise \code{"flow"}
#'           (fluorescence cytometry). The result is stored in
#'           \code{obj$collection_instrument}, and an update timestamp is
#'           stored in \code{obj$object_history}.
#'   }
#'
#' @param fcs_join_obj
#'   An FCSimple analysis object (list) returned by \code{FCSimple::fcs_join()}.
#'   Must contain the fields \code{data}, \code{source}, and \code{run_date}.
#'
#' @details
#'   \strong{Metadata preservation:} Existing metadata is kept intact if all
#'   four conditions hold: (1) \code{metadata} is a \code{data.frame},
#'   (2) it contains a \code{patient_ID} column, (3) every value in
#'   \code{patient_ID} appears in \code{source}, (4) every unique value in
#'   \code{source} appears in \code{patient_ID}, and (5) the metadata has
#'   exactly one row per unique \code{patient_ID}. If \code{run_date} is
#'   missing from the otherwise-valid metadata it is joined from the
#'   cell-level \code{source}/\code{run_date} vectors. If any condition
#'   fails, metadata is rebuilt from scratch.
#'
#'   \strong{Version capture:} Python availability is tested with
#'   \code{system("python --version")}. If Python is found, the installed
#'   package list is retrieved via \code{pip list} (or \code{pip3 list} as a
#'   fallback) and stored as a two-column \code{data.frame} (\code{Package},
#'   \code{Version}). If Python is not available, both \code{Python} and
#'   \code{pip_list} are set to \code{"none"}.
#'
#'   \strong{Instrument detection:} Because CyTOF/mass-cytometry ion counts are
#'   non-negative whilst flow-cytometry fluorescence values can be negative
#'   (after compensation), the sign of \code{min(obj$data)} is used as a
#'   heuristic: \eqn{\geq 0} → \code{"cytof"}, \eqn{< 0} → \code{"flow"}.
#'
#' @return
#'   The input \code{fcs_join_obj} list with the following fields added or
#'   updated:
#'   \describe{
#'     \item{\code{versions}}{A list with elements \code{R} (R version),
#'       \code{Python} (Python version string or \code{"none"}),
#'       \code{Rsession} (output of \code{sessionInfo()}), and
#'       \code{pip_list} (pip package \code{data.frame} or \code{"none"}).}
#'     \item{\code{metadata}}{A \code{data.frame} with at least columns
#'       \code{patient_ID} and \code{run_date}, one row per unique sample.
#'       Preserved from the input if it already satisfies the coverage and
#'       uniqueness requirements; rebuilt from \code{source} and
#'       \code{run_date} otherwise.}
#'     \item{\code{collection_instrument}}{Character string: \code{"cytof"} or
#'       \code{"flow"}.}
#'     \item{\code{object_history}}{Character string recording the update
#'       timestamp, e.g. \code{"updated: 2024-06-01 14:32:07"}.}
#'   }
#'
#' @examples
#' \dontrun{
#'   files <- list.files("data/fcs", "\\.fcs$", full.names = TRUE)
#'   joined <- FCSimple::fcs_join(files)
#'
#'   # Update: patient_ID is taken directly from the source field
#'   updated <- FCSimple::fcs_update(joined)
#'
#'   # Inspect added fields
#'   updated$collection_instrument  # "flow" or "cytof"
#'   updated$object_history         # "updated: <timestamp>"
#'   head(updated$metadata)         # patient_ID, run_date
#'   updated$versions$R             # R version used
#' }
#'
#' @seealso
#'   \code{\link[FCSimple]{fcs_join}}, \code{\link[FCSimple]{fcs_gating_object}}
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
  if(mean(c("data", "source", "run_date") %in% names(fcs_join_obj))!=1) {
    stop('fcs_join_obj must contain: "data", "source", and "run_date"')
  }
  unique_sources <- unique(fcs_join_obj[['source']])
  existing_meta  <- fcs_join_obj[['metadata']]
  meta_is_valid  <- !is.null(existing_meta) &&
                    is.data.frame(existing_meta) &&
                    'patient_ID' %in% names(existing_meta) &&
                    all(existing_meta[['patient_ID']] %in% unique_sources) &&
                    all(unique_sources %in% existing_meta[['patient_ID']]) &&
                    length(unique(existing_meta[['patient_ID']])) == length(unique_sources) &&
                    nrow(existing_meta) == length(unique(existing_meta[['patient_ID']]))
  if(meta_is_valid) {
    # Metadata already fully covers all sources; preserve it, add run_date if absent
    if(!'run_date' %in% names(fcs_join_obj[['metadata']])) {
      rd_map <- unique(data.frame(patient_ID = fcs_join_obj[['source']],
                                  run_date   = fcs_join_obj[['run_date']],
                                  stringsAsFactors = FALSE))
      rd_map <- rd_map[!duplicated(rd_map$patient_ID), ]
      fcs_join_obj[['metadata']][['run_date']] <-
        rd_map[['run_date']][match(fcs_join_obj[['metadata']][['patient_ID']], rd_map[['patient_ID']])]
    }
  } else {
    # Build metadata from scratch
    base_metadata <- data.frame(patient_ID = fcs_join_obj[['source']],
                                run_date   = fcs_join_obj[['run_date']],
                                stringsAsFactors = FALSE)
    fcs_join_obj[['metadata']] <- base_metadata[!duplicated(base_metadata$patient_ID), ]
  }
  
  # Auto-detect instrument type based on data
  detected_instrument <- if(min(fcs_join_obj$data, na.rm = TRUE) >= 0) "cytof" else "flow"
  
  fcs_join_obj[['collection_instrument']] <- detected_instrument
  fcs_join_obj[['object_history']] <- paste0("updated: ",Sys.time())
  return(fcs_join_obj)
}
