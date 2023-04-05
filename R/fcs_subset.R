fcs_subset <- function(fcs_join_obj,
                       subset_by = c("cluster","source"),
                       subset_cluster_algorithm = NA,
                       subset_values)
{
  if(subset_by=="source") {
    get_index <- which(fcs_join_obj[["source"]] %in% subset_values)
    if(length(get_index)==0) {
      stop(paste0("Value requested not found in source. Use one or more of: ",
                  paste0(names(table(fcs_join_obj[["source"]])), collapse = ", ")))
    }
    fcs_new_obj <- list(data = fcs_join_obj[["data"]][get_index,],
                        source = fcs_join_obj[["source"]][get_index])
    if("run_date" %in% names(fcs_join_obj)) {
      fcs_new_obj[["run_date"]] <- fcs_join_obj[["run_date"]][get_index]
    }
  } else if(subset_by=="cluster") {
    cluster_nums <- fcs_join_obj[[tolower(subset_cluster_algorithm)]][["clusters"]]
    get_index <- which(cluster_nums %in% subset_values)
    fcs_new_obj <- list(data = fcs_join_obj[["data"]][get_index,],
                        source = fcs_join_obj[["source"]][get_index])
    if("run_date" %in% names(fcs_join_obj)) {
      fcs_new_obj[["run_date"]] <- fcs_join_obj[["run_date"]][get_index]
    }
  } else {
    stop("error in argument 'subset_by': only subettable by cluster or source")
  }
  fcs_new_obj[["subset_on"]] <- list(subset_by = subset_by,
                                     subset_cluster_algorithm = subset_cluster_algorithm,
                                     subset_values = subset_values,
                                     source_object_indices = get_index)
  return(fcs_new_obj)
}

fcs_remove_parameters <- function(fcs_join_obj,
                                  remove_parameters)
{
  rm_param <- which(colnames(fcs_join_obj[["data"]]) %in% remove_parameters)
  if(length(rm_param)==0) {
    stop("error in argument 'remove_parameters': none of the requested parameters were found")
  } else {
    new_obj <- list(data = fcs_join_obj[["data"]][,-rm_param],
                    source = fcs_join_obj[["source"]])
    if("run_date" %in% names(fcs_join_obj)) {
      new_obj[["run_date"]] <- fcs_join_obj[["run_date"]]
    }
    return(new_obj)
  }
}
