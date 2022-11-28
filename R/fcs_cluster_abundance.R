fcs_calculate_abundance <- function(fcs_join_obj,
                                    report_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                    report_as = c("frequency", "fraction"))
{
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
      fval <- mean(tmp_numbers==as.numeric(colnames(frequency_matrix)[j]))
      if(report_as=="frequency") {
        frequency_matrix[i,j] <- fval * 100
      } else if(report_as=="fraction") {
        frequency_matrix[i,j] <- fval
      }
    }
  }
  fcs_join_obj[[report_algorithm]][["abundance"]] <- frequency_matrix
  return(fcs_join_obj)
}

fcs_report_abundance <- function(fcs_join_obj,
                                 report_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                 outdir = getwd())
{
  abundance_values <- fcs_join_obj[[tolower(report_algorithm)]][["abundance"]]
  row.names(abundance_values) <- gsub("^.+/|.fcs$","",row.names(abundance_values))
  outdir <- gsub("/$","",outdir)
  write.csv(x = abundance_values, file = paste0(outdir,"/",report_algorithm,"_cluster_abundance_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S")), row.names = TRUE)
  return(abundance_values)
}
