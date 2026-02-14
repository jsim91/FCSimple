#' @title Add Custom Metadata to an FCSimple Object
#'
#' @description
#'   Merges user-supplied metadata into the `metadata` slot of an FCSimple
#'   analysis object created by `FCSimple::fcs_join()`. The custom metadata
#'   is joined on the `patient_ID` column using `merge()`, allowing you to
#'   add clinical, demographic, or other external information to your flow
#'   or mass cytometry dataset.
#'
#' @param fcs_join_obj
#'   An FCSimple analysis object (list) returned by `FCSimple::fcs_join()`.
#'   Must contain a `metadata` element with a `patient_ID` column.
#'
#' @param custom_metadata
#'   A data frame containing additional metadata to merge. Must have a
#'   `patient_ID` column that matches values in `fcs_join_obj$metadata$patient_ID`.
#'   Additional columns will be added to the object's metadata.
#'
#' @param suffix
#'   Character vector of length 2; suffixes to append to duplicate column
#'   names (excluding `patient_ID`). Default `c(".obj", ".custom")`.
#'   Passed to `merge()`.
#'
#' @details
#'   The function validates that:
#'   - `fcs_join_obj` is a list with a `metadata` element
#'   - Both `fcs_join_obj$metadata` and `custom_metadata` have a `patient_ID` column
#'   - There is at least one matching `patient_ID` between the two datasets
#'
#'   The merge is performed using `base::merge()`, which performs a SQL-like join.
#'   After merging, the updated metadata is stored back in `fcs_join_obj$metadata`,
#'   and a timestamped entry is appended to `object_history`.
#'
#' @return
#'   The input `fcs_join_obj` with an updated `metadata` element containing
#'   the merged data, and an updated `object_history` recording the operation.
#'
#' @examples
#' \dontrun{
#'   # Load and join FCS files
#'   files <- list.files("data/fcs", "\\.fcs$", full.names = TRUE)
#'   joined <- FCSimple::fcs_join(files)
#'
#'   # Create custom metadata with patient_ID and clinical info
#'   clinical_data <- data.frame(
#'     patient_ID = c("patient1", "patient2", "patient3"),
#'     age = c(45, 52, 38),
#'     treatment = c("A", "B", "A"),
#'     response = c("CR", "PR", "NR")
#'   )
#'
#'   # Add custom metadata to the object
#'   joined <- FCSimple::fcs_add_metadata(joined, clinical_data)
#'
#'   # View the updated metadata
#'   head(joined$metadata)
#' }
#'
#' @seealso
#'   FCSimple::fcs_join, FCSimple::fcs_update, base::merge
#'
#' @export
fcs_add_metadata <- function(fcs_join_obj,
                             custom_metadata,
                             suffix = c(".obj", ".custom")) {
  
  # Validate inputs
  if (!is.list(fcs_join_obj)) {
    stop("fcs_join_obj must be a list (FCSimple analysis object)")
  }
  
  if (!"metadata" %in% names(fcs_join_obj)) {
    stop("fcs_join_obj must contain a 'metadata' element. Use FCSimple::fcs_join() to create a proper object.")
  }
  
  if (!is.data.frame(custom_metadata)) {
    stop("custom_metadata must be a data frame")
  }
  
  # Check for patient_ID column in both datasets
  if (!"patient_ID" %in% colnames(fcs_join_obj$metadata)) {
    stop("fcs_join_obj$metadata must have a 'patient_ID' column")
  }
  
  if (!"patient_ID" %in% colnames(custom_metadata)) {
    stop("custom_metadata must have a 'patient_ID' column for merging")
  }
  
  # Check for at least one matching patient_ID
  matching_ids <- intersect(
    fcs_join_obj$metadata$patient_ID,
    custom_metadata$patient_ID
  )
  
  if (length(matching_ids) == 0) {
    stop("No matching patient_ID values found between fcs_join_obj$metadata and custom_metadata")
  }
  
  # Inform user about the merge
  n_matches <- length(matching_ids)
  n_total_obj <- nrow(fcs_join_obj$metadata)
  n_total_custom <- nrow(custom_metadata)
  
  message(sprintf(
    "Merging metadata: %d of %d patient_ID values in fcs_join_obj match %d of %d in custom_metadata",
    n_matches, n_total_obj, n_matches, n_total_custom
  ))
  
  # Perform the merge
  merged_metadata <- merge(
    x = fcs_join_obj$metadata,
    y = custom_metadata,
    by = "patient_ID",
    all.x = TRUE,
    all.y = FALSE,
    suffix = suffix
  )
  
  # Update the object
  fcs_join_obj$metadata <- merged_metadata
  
  # Update object history
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  history_entry <- paste0("fcs_add_metadata: ", timestamp, " - Added ", 
                          ncol(custom_metadata) - 1, " custom metadata column(s)")
  
  if ("object_history" %in% names(fcs_join_obj)) {
    fcs_join_obj$object_history <- c(fcs_join_obj$object_history, history_entry)
  } else {
    fcs_join_obj$object_history <- history_entry
  }
  
  return(fcs_join_obj)
}
