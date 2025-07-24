#' @title Calculate Cluster Abundance Across Samples
#'
#' @description
#'   Computes per‐sample cluster abundances for a previously clustered flow
#'   cytometry object. Stores a matrix of frequencies (0–100) or fractions
#'   (0–1) under the chosen algorithm’s `abundance` element, and appends an
#'   entry to `object_history`.
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
#'   `"leiden"`, `"flowsom"`, `"louvain"`, `"phenograph"`, or `"git"`.
#'
#' @param report_as
#'   Character; type of abundance metric. `"frequency"` (0–100) or
#'   `"fraction"` (0–1). Defaults to `"frequency"`.
#'
#' @details
#'   The function will:
#'   1. Validate that the specified algorithm and `source` exist.
#'   2. Tabulate the proportion of cells in each cluster, per sample.
#'   3. Multiply by 100 if `report_as = "frequency"`.
#'   4. Store the result in `fcs_join_obj[[report_algorithm]][["abundance"]]`.
#'   5. Append a timestamped note to `object_history`.
#'
#' @return
#'   The input `fcs_join_obj`, augmented with:
#'   - `<algorithm>$abundance`: a numeric matrix (samples × clusters)
#'   - updated `object_history` entry
#'
#' @examples
#' \dontrun{
#'   joined <- FCSimple::fcs_join(list(ff1, ff2))
#'   clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#'
#'   # Get percentage abundance for Leiden clusters
#'   out <- FCSimple::fcs_calculate_abundance(
#'     clustered,
#'     report_algorithm = "leiden",
#'     report_as = "frequency"
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_cluster, FCSimple::fcs_report_abundance
#'
#' @importFrom utils write.csv
#' @export
fcs_calculate_abundance <- function(fcs_join_obj,
                                    report_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                    report_as = c("frequency", "fraction"))
{
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  }
  if(!report_algorithm %in% names(fcs_join_obj)) {
    stop("error in names of fcs_join_obj: has data been clustered using a supported algorithm. See ?cluster.")
  }
  if(!"source" %in% names(fcs_join_obj)) {
    stop("error in names of fcs_join_obj: could not trace cells back to their origin. See ?fcs_join.")
  }
  if(length(report_as)>2) {
    report_as <- tolower(report_as[1])
    if(!report_as %in% c("frequency", "fraction")) {
      stop("error in argument 'report_as': use either 'frequency' (range 0-100) or 'fraction' (range 0-1)")
    }
  }
  cluster_numbers <- fcs_join_obj[[which(tolower(names(fcs_join_obj))==tolower(report_algorithm))]][[1]]
  cluster_source <- fcs_join_obj[["source"]]
  usrc <- unique(cluster_source); uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  frequency_matrix <- matrix(data = NA, nrow = length(usrc), ncol = length(uclus))
  row.names(frequency_matrix) <- usrc; colnames(frequency_matrix) <- uclus
  for(i in 1:nrow(frequency_matrix)) {
    tmp_numbers <- cluster_numbers[which(cluster_source==row.names(frequency_matrix)[i])]
    for(j in 1:ncol(frequency_matrix)) {
      fval <- mean(tmp_numbers==colnames(frequency_matrix)[j])
      if(report_as=="frequency") {
        frequency_matrix[i,j] <- fval * 100
      } else if(report_as=="fraction") {
        frequency_matrix[i,j] <- fval
      }
    }
  }
  fcs_join_obj[[report_algorithm]][["abundance"]] <- frequency_matrix
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  } else {
    fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(tolower(report_algorithm)," abundance calculated: ",Sys.time()))
  }
  return(fcs_join_obj)
}

#' @title Report Cluster Abundance to CSV File
#'
#' @description
#'   Exports the cluster‐abundance matrix (samples × clusters) to a CSV in the
#'   specified directory, and returns the matrix invisibly.
#'
#' @param fcs_join_obj
#'   A list with a `[[algorithm]][["abundance"]]` matrix, as produced by
#'   FCSimple::fcs_calculate_abundance().
#'
#' @param report_algorithm
#'   Character; name of the abundance matrix to export. One of
#'   `"leiden"`, `"flowsom"`, `"louvain"`, `"phenograph"`, or `"git"`.
#'
#' @param outdir
#'   Character; file path to an existing directory. Defaults to `getwd()`.
#'
#' @return
#'   Invisibly returns the abundance matrix (numeric matrix with sample rows
#'   and cluster‐ID columns).
#'
#' @examples
#' \dontrun{
#'   # Assume 'clustered' has abundance computed
#'   abundance_mat <- FCSimple::fcs_calculate_abundance(clustered)
#'
#'   # Write to your working directory
#'   FCSimple::fcs_report_abundance(
#'     clustered,
#'     report_algorithm = "leiden",
#'     outdir = "~/my_analysis/results"
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_calculate_abundance
#'
#' @importFrom utils write.csv
#' @export
fcs_report_abundance <- function(fcs_join_obj,
                                 report_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                 outdir = getwd())
{
  abundance_values <- fcs_join_obj[[tolower(report_algorithm)]][["abundance"]]
  row.names(abundance_values) <- gsub("^.+/|.fcs$","",row.names(abundance_values))
  outdir <- gsub("/$","",outdir)
  write.csv(x = abundance_values, file = paste0(outdir,"/",report_algorithm,"_cluster_abundance_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".csv"), row.names = TRUE)
  return(abundance_values)
}
