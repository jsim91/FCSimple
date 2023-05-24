fcs_write.FCS <- function(fcs_join_obj,
                          fcs_name = "fcs_out",
                          data_format = c("raw","transformed"),
                          include_reductions = c("UMAP","tSNE"),
                          include_clusterings = c("leiden","flowsom","louvain","phenograph","git"),
                          outdir = getwd(),
                          include_timestamp = TRUE)
{
  require(flowCore)

  if(length(data_format)==2) {
    data_format <- "raw"
  }
  if(data_format=="raw") {
    data_incl <- fcs_join_obj[["raw"]]
  } else if(data_format=="transformed") {
    data_incl <- fcs_join_obj[["data"]]
  } else {
    stop("error in argument 'data_format': use either 'raw' or 'transformed'")
  }

  fcs_include <- vector("list", length = 3L)
  fcs_include[[1]] <- data_incl

  if(length(include_reductions)!=0) {
    hold_reductions <- vector("list", length = length(include_reductions))
    for(i in 1:length(include_reductions)) {
      if(tolower(include_reductions[i]) %in% names(fcs_join_obj)) {
        hold_reductions[[i]] <- as.matrix(fcs_join_obj[[tolower(include_reductions[i])]][["coordinates"]])
      } else {
        warning(paste0("results of reduction ",tolower(include_reductions[i])," not found in object and will not be included"))
        hold_reductions[[i]] <- NULL
      }
    }
    rm_reduction <- c()
    for(i in 1:length(hold_reductions)) {
      if(is.null(hold_reductions[[i]])) {
        rm_reduction <- append(rm_reduction, i)
      }
    }
    if(length(rm_reduction)!=0) {
      if(length(rm_reduction)==1) {
        pared_reduction <- hold_reductions[[-rm_reduction]]
      } else if(length(rm_reduction)>1) {
        pared_reduction <- hold_reductions[-rm_reduction]
      }
    } else if(length(rm_reduction)==0) {
      pared_reduction <- hold_reductions
    }
    for(i in 1:length(pared_reduction)) {
      if(i==1) {
        add_reduction <- pared_reduction[[i]]
      } else {
        add_reduction <- cbind(add_reduction, pared_reduction[[i]])
      }
    }
    fcs_include[[2]] <- as.matrix(add_reduction)
  } else {
    fcs_include[[2]] <- NULL
  }
  if(length(include_clusterings)!=0) {
    cluster_numbers <- vector("list", length = length(include_clusterings))
    names(cluster_numbers) <- paste0(tolower(include_clusterings),"_cluster")
    for(i in 1:length(include_clusterings)) {
      if(tolower(include_clusterings[i]) %in% names(fcs_join_obj)) {
        cluster_numbers[[i]] <- fcs_join_obj[[tolower(include_clusterings[i])]][["clusters"]]
      } else {
        warning(paste0("results of clustering ",tolower(include_clusterings[i])," not found in object and will not be included"))
        cluster_numbers[[i]] <- NULL
      }
    }
    rm_cluster <- c()
    for(i in 1:length(cluster_numbers)) {
      if(is.null(cluster_numbers[[i]])) {
        rm_reduction <- append(rm_reduction, i)
      }
    }
    if(length(rm_reduction)!=0) {
      if(length(rm_reduction)==1) {
        pared_cluster <- cluster_numbers[[-rm_reduction]]
      } else if(length(rm_reduction)>1) {
        pared_cluster <- cluster_numbers[-rm_reduction]
      }
    } else if(length(rm_reduction)==0) {
      pared_cluster <- cluster_numbers
    }
    for(i in 1:length(pared_cluster)) {
      if(i==1) {
        tmp_cluster_mat <- data.frame(var1 = pared_cluster[[i]])
        colnames(tmp_cluster_mat)[1] <- paste0(names(pared_clusters)[i],"_cluster")
        add_cluster <- tmp_cluster_mat
      } else {
        tmp_cluster_mat <- data.frame(var1 = pared_cluster[[i]])
        colnames(tmp_cluster_mat)[1] <- paste0(names(pared_clusters)[i],"_cluster")
        add_cluster <- cbind(add_cluster, tmp_cluster_mat)
      }
    }
    fcs_include[[3]] <- as.matrix(add_cluster)
  } else {
    fcs_include[[3]] <- NULL
  }
  which_populated <- as.numeric(which(sapply(fcs_include, function(arg1) return(!is.null(arg1)))))
  if(length(which_populated)==0) {
    stop("unable to return any values.. verify that the input object has a 'data' entry")
  } else if(length(which_populated)==1) {
    final_return <- fcs_include[[which_populated]]
  } else if(length(which_populated)>1) {
    for(i in 1:length(which_populated)) {
      if(i==1) {
        final_return <- fcs_include[[which_populated[[i]]]]
      } else {
        final_return <- cbind(final_return, fcs_include[[which_populated[[i]]]])
      }
    }
  }
  final_return <- as.matrix(final_return) # structure should be a matrix, but as.matrix is used just in case
  out_ff <- new("flowFrame", exprs = final_return)
  if(include_timestamp) {
    write.FCS(x = out_ff, filename = paste0(gsub("\\/$","",outdir),"/",gsub("\\.fcs$","",fcs_name),"_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".fcs"))
  } else {
    write.FCS(x = out_ff, filename = paste0(gsub("\\/$","",outdir),"/",gsub("\\.fcs$","",fcs_name),".fcs"))
  }
}
