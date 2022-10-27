fcs_plot_distributions <- function(fcs_join_obj,
                                   separate_by = c("none", "date"),
                                   plot_element = c("cluster","total"),
                                   plot_algorithm = c("leiden","flowsom","louvain","phenograph"))
{
  if(separate_by=="date") {
    if(!"run_date" %in% names(fcs_join_obj)){
      print("Unable to find run date. Using separate_by = 'none' instead.")
      obj_data <- fcs_join_obj[["data"]]
    } else {
      obj_data <- fcs_join_obj[["data"]]
      batch <- fcs_join_obj[["run_date"]]
      obj_data$date <- batch
    }
  }
  if("batch" %in% obj_data) {
    data_split <- vector("list", length = ncol(obj_data)-1)
    names(data_split) <- colnames(obj_data)[1:(ncol(obj_data)-1)]
    for(i in 1:length(data_split)) {
      data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])],
                                    batch = obj_data[,"batch"])
    }
  } else {
    data_split <- vector("list", length = ncol(obj_data))
    data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])])
  }
  if(tolower(plot_element) == "cluster") {
    if(!tolower(plot_algorithm) %in% c("leiden","flowsom","louvain","phenograph")) {
      stop(paste0("error in argument 'plot_algorithm': clusters not found for ",tolower(plot_algorithm)," not found in 'fcs_join_obj'"))
    }
    for(i in 1:length(data_split)) {
      data_split[[i]]$cluster <- fcs_join_obj[[tolower(plot_algorithm)]]$clusters
    }
  }
}

fcs_plot_help <- function(input_data)
{

}
