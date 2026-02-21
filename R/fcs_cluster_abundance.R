#' @title Calculate Cluster Abundance Across Samples
#'
#' @description
#'   Computes per‐sample cluster abundances for a previously clustered flow
#'   cytometry object. Stores the result under the chosen algorithm's
#'   `frequency`, `fraction`, or `counts` element (matching `report_as`),
#'   and appends a timestamped entry to `object_history`.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join() and FCSimple::fcs_cluster(),
#'   containing at minimum:
#'   - `source`: a character vector of original sample identifiers
#'   - `<algorithm>`: a list whose first element is a vector of cluster IDs
#'   - `object_history`: an optional history log
#'
#' @param report_algorithm
#'   Character; name of the clustering result to summarize. One of
#'   `"leiden"`, `"flowsom"`, `"louvain"`, or `"phenograph"`.
#'
#' @param report_as
#'   Character; type of abundance metric. One of
#'   `"frequency"` (0–100), `"fraction"` (0–1), or `"count"` (raw cell counts).
#'   Defaults to `"frequency"`.
#'
#' @param return_abundance
#'   Logical; if `TRUE`, the function returns the abundance matrix directly
#'   and does not modify `fcs_join_obj`. Defaults to `FALSE`.
#'
#' @details
#'   The function will:
#'   1. Validate that the specified algorithm and `source` exist.
#'   2. If `report_as = "count"`, tabulate raw cell counts per cluster, per sample.
#'      Otherwise, compute proportions of cells in each cluster per sample.
#'   3. If `report_as = "frequency"`, multiply proportions by 100.
#'   4. If `return_abundance = TRUE`, return the matrix directly.
#'      Otherwise, store it in `fcs_join_obj[[report_algorithm]][[report_as]]`
#'      (i.e. `[["frequency"]]`, `[["fraction"]]`, or `[["counts"]]`).
#'   5. If not returning directly, append a timestamped note to
#'      `object_history`.
#'
#' @return
#'   If `return_abundance = TRUE`, a numeric matrix (samples × clusters)
#'   containing counts, frequencies, or fractions. Otherwise, the input
#'   `fcs_join_obj`, augmented with:
#'   - `<algorithm>$frequency`: a numeric matrix (samples × clusters) when `report_as = "frequency"`
#'   - `<algorithm>$fraction`: a numeric matrix (samples × clusters) when `report_as = "fraction"`
#'   - `<algorithm>$counts`: a numeric matrix (samples × clusters) when `report_as = "count"`
#'   - updated `object_history`
#'
#' @examples
#' \dontrun{
#'   joined <- FCSimple::fcs_join(list(ff1, ff2))
#'   clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#'
#'   # Get percentage abundance for Leiden clusters
#'   out1 <- FCSimple::fcs_calculate_abundance(
#'     clustered,
#'     report_algorithm = "leiden",
#'     report_as = "frequency"
#'   )
#'
#'   # Get raw cell counts for Leiden clusters, return matrix directly
#'   counts_mat <- FCSimple::fcs_calculate_abundance(
#'     clustered,
#'     report_algorithm = "leiden",
#'     report_as = "count",
#'     return_abundance = TRUE
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_cluster, FCSimple::fcs_report_abundance
#'
#' @importFrom utils write.csv
#' @export
fcs_calculate_abundance <- function(fcs_join_obj,
                                    report_algorithm = c("leiden","flowsom","louvain","phenograph"),
                                    report_as = c("frequency", "fraction", "count"), return_abundance = FALSE)
{
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_audit() on the object.")
  }
  if(!report_algorithm %in% names(fcs_join_obj)) {
    stop("error in names of fcs_join_obj: has data been clustered using a supported algorithm. See ?cluster.")
  }
  if(!"source" %in% names(fcs_join_obj)) {
    stop("error in names of fcs_join_obj: could not trace cells back to their origin. See ?fcs_join.")
  }
  if(length(report_as)>2) {
    report_as <- tolower(report_as[1])
    warning(paste0("more than one 'report_as' specified... continuing with ",report_as[1]))
    if(!report_as %in% c("frequency", "fraction","count")) {
      stop("error in argument 'report_as': use either 'frequency' (range 0-100), 'fraction' (range 0-1), or 'count' (counts table)")
    }
  }
  cluster_numbers <- fcs_join_obj[[which(tolower(names(fcs_join_obj))==tolower(report_algorithm))]][[1]]
  cluster_source <- fcs_join_obj[["source"]]
  if(report_as=='count') {
    abundance_cross <- table(cluster_source, cluster_numbers)
    abundance_matrix <- as.matrix(as.data.frame.matrix(abundance_cross))
    if(return_abundance) {
      return(abundance_matrix)
    } else {
      fcs_join_obj[[report_algorithm]][["counts"]] <- abundance_matrix
    }
  } else {
    usrc <- unique(cluster_source); uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]
    abundance_matrix <- matrix(data = NA, nrow = length(usrc), ncol = length(uclus))
    row.names(abundance_matrix) <- usrc; colnames(abundance_matrix) <- uclus
    for(i in 1:nrow(abundance_matrix)) {
      tmp_numbers <- cluster_numbers[which(cluster_source==row.names(abundance_matrix)[i])]
      for(j in 1:ncol(abundance_matrix)) {
        fval <- mean(tmp_numbers==colnames(abundance_matrix)[j])
        if(report_as=="frequency") {
          abundance_matrix[i,j] <- fval * 100
        } else if(report_as=="fraction") {
          abundance_matrix[i,j] <- fval
        }
      }
    }
    if(return_abundance) {
      return(abundance_matrix)
    } else {
      fcs_join_obj[[report_algorithm]][[report_as]] <- abundance_matrix
    }
  }
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_audit() on the object.")
  } else {
    fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(tolower(report_algorithm)," ",report_as," calculated: ",Sys.time()))
  }
  return(fcs_join_obj)
}

#' @title Report Cluster Frequency or Fraction to CSV File
#'
#' @description
#' Exports a cluster frequency or fraction matrix (samples × clusters) to a
#' CSV file in the specified directory and returns the matrix invisibly. By
#' default, the function prepends one timestamp and can optionally add a
#' second timestamp or custom string to the filename.
#'
#' @param fcs_join_obj A list containing an `[[algorithm]][[report_as]]`
#'   numeric matrix (i.e. `[["frequency"]]` or `[["fraction"]]`), as produced
#'   by `FCSimple::fcs_calculate_abundance()`.
#'
#' @param report_algorithm Character scalar. Name of the clustering result to
#'   export; must be one of `"leiden"`, `"flowsom"`, `"louvain"`,
#'   or `"phenograph"`. Defaults to `"leiden"`.
#'
#' @param report_as Character scalar. Which stored matrix to export: either
#'   `"frequency"` (0–100, the default) or `"fraction"` (0–1). Must match
#'   what was calculated by `fcs_calculate_abundance()`.
#'
#' @param outdir Character scalar. Path to an existing directory where the
#'   CSV will be written. Trailing slashes are removed internally.
#'   Defaults to `getwd()`.
#'
#' @param add_timestamp Logical scalar. If `TRUE` (the default), a second
#'   timestamp (`%Y-%m-%d_%H%M%S`) is appended to the filename after the
#'   initial timestamp.
#'
#' @param append_file_string Character scalar or `NA`. If non‐`NA`, this
#'   string is appended (after any timestamps) to the filename. Defaults to
#'   `NA`, which means no extra string is added.
#'
#' @details
#' The output filename is built in three parts:
#'     1. `<report_algorithm>_cluster_frequency_<timestamp1>` (or `fraction`)
#'     2. `_<timestamp2>` (only if `add_timestamp = TRUE`)
#'     3. `_<append_file_string>` (only if `append_file_string` is not `NA`)
#'
#' All timestamps use the format `%Y-%m-%d_%H%M%S`. The final suffix `.csv`
#' is then added.
#'
#' @return
#' Invisibly returns the same frequency or fraction matrix that was written
#' to disk (a numeric matrix with sample rows and cluster‐ID columns).
#'
#' @examples
#' \dontrun{
#' # Calculate cluster frequencies first
#' clustered <- FCSimple::fcs_calculate_abundance(clustered, report_algorithm = "leiden", report_as = "frequency")
#'
#' # Write frequency matrix to disk
#' FCSimple::fcs_report_abundance(
#'   clustered,
#'   report_algorithm   = "leiden",
#'   report_as          = "frequency",
#'   outdir             = "~/my_analysis/results",
#'   add_timestamp      = FALSE,
#'   append_file_string = "v1"
#' )
#' }
#'
#' @seealso
#' FCSimple::fcs_calculate_abundance
#'
#' @importFrom utils write.csv
#' @export
fcs_report_abundance <- function(fcs_join_obj,
                                 report_algorithm = c("leiden","flowsom","louvain","phenograph"),
                                 report_as = c("frequency", "fraction"),
                                 outdir = getwd(), add_timestamp = TRUE, append_file_string = NA)
{
  report_as <- tolower(report_as[1])
  abundance_values <- fcs_join_obj[[tolower(report_algorithm)]][[report_as]]
  row.names(abundance_values) <- gsub("^.+/|.fcs$","",row.names(abundance_values))
  outdir <- gsub("/$","",outdir)
  fname <- paste0(report_algorithm,"_cluster_",report_as,"_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
  if(add_timestamp) {
    fname <- paste0(fname,'_',strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
  }
  if(!is.na(append_file_string)) {
    fname <- paste0(fname,'_',append_file_string)
  }
  fname <- paste0(fname, '.csv')
  write.csv(x = abundance_values, file = file.path(outdir,fname), row.names = TRUE)
}
