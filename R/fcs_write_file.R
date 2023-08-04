fcs_write.FCS <- function(fcs_join_obj,
                          fcs_name = "fcs_out",
                          data_format = c("raw","transformed"),
                          include_reductions = c("UMAP","tSNE"),
                          include_clusterings = c("leiden","flowsom","louvain","phenograph","git"),
                          subset_rows = "all", # either "all" or a numeric vector of positions to write to file
                          outdir = getwd(),
                          include_timestamp = TRUE)
{
  require(flowCore)

  if(length(data_format)==2) {
    data_format <- "raw"
  }
  if(data_format=="raw") {
    if(!"raw" %in% names(fcs_join_obj)) {
      warning("'raw' format not found in object. Consider adding a 'raw' element to the list object with reverse transformed values. Using transformed values stored in $data element.")
      data_incl <- fcs_join_obj[["data"]]
    } else {
      data_incl <- fcs_join_obj[["raw"]]
    }
  } else if(data_format=="transformed") {
    data_incl <- fcs_join_obj[["data"]]
  } else {
    stop("error in argument 'data_format': use either 'raw' or 'transformed'")
  }

  fcs_include <- list('1' = 1, '2' = 2, '3' = 3)
  fcs_include[[1]] <- data_incl

  if(length(include_reductions)!=0) {
    hold_reductions <- vector("list", length = length(include_reductions))
    for(i in 1:length(include_reductions)) {
      if(tolower(include_reductions[i]) %in% names(fcs_join_obj)) {
        hold_reductions[[i]] <- as.matrix(fcs_join_obj[[tolower(include_reductions[i])]][["coordinates"]])
      } else {
        warning(paste0("results of reduction ",tolower(include_reductions[i])," not found in object and will not be included"))
        hold_reductions[[i]] <- NA
      }
    }
    rm_reduction <- c()
    for(i in 1:length(hold_reductions)) {
      if(class(hold_reductions[[i]])[1]=="logical") {
        if(is.na(hold_reductions[[i]][1]))
          rm_reduction <- append(rm_reduction, i)
      }
    }
    pared_reduction <- hold_reductions
    if(length(rm_reduction)!=0) {
      if(length(rm_reduction)==1) {
        pared_reduction[[rm_reduction]] <- NULL
      } else if(length(rm_reduction)>1) {
        pared_reduction[rm_reduction] <- NULL
      }
    } else if(length(rm_reduction)==0) {
      pared_reduction <- hold_reductions
    }
    if(length(pared_reduction)==0) {
      fcs_include[[2]] <- NA
    } else {
      for(i in 1:length(pared_reduction)) {
        if(i==1) {
          add_reduction <- pared_reduction[[i]]
        } else {
          add_reduction <- cbind(add_reduction, pared_reduction[[i]])
        }
      }
      fcs_include[[2]] <- as.matrix(add_reduction)
    }
  } else {
    fcs_include[[2]] <- NA
  }
  if(length(include_clusterings)!=0) {
    cluster_numbers <- vector("list", length = length(include_clusterings))
    names(cluster_numbers) <- paste0(tolower(include_clusterings),"_cluster")
    for(i in 1:length(include_clusterings)) {
      if(tolower(include_clusterings[i]) %in% names(fcs_join_obj)) {
        cluster_numbers[[i]] <- fcs_join_obj[[tolower(include_clusterings[i])]][["clusters"]]
      } else {
        warning(paste0("results of clustering ",tolower(include_clusterings[i])," not found in object and will not be included"))
        cluster_numbers[[i]] <- NA
      }
    }
    rm_cluster <- c()
    for(i in 1:length(cluster_numbers)) {
      if(class(cluster_numbers[[i]])[1]=="logical") {
        if(is.na(cluster_numbers[[i]][1]))
          rm_cluster <- append(rm_cluster, i)
      }
    }
    pared_cluster <- cluster_numbers
    if(length(rm_cluster)!=0) {
      if(length(rm_cluster)==1) {
        pared_cluster[[rm_cluster]] <- NULL
      } else if(length(rm_cluster)>1) {
        pared_cluster[rm_cluster] <- NULL
      }
    } else if(length(rm_cluster)==0) {
      pared_cluster <- cluster_numbers
    }
    if(length(pared_cluster)==0) {
      fcs_include[[3]] <- NA
    } else {
      for(i in 1:length(pared_cluster)) {
        if(i==1) {
          tmp_cluster <- data.frame(var1 = pared_cluster[[i]])
          colnames(tmp_cluster)[i] <- names(pared_cluster)[i]
          add_cluster <- as.matrix(tmp_cluster)
        } else {
          tmp_cluster <- data.frame(var1 = pared_cluster[[i]])
          colnames(tmp_cluster)[i] <- names(pared_cluster)[i]
          add_cluster <- cbind(add_cluster, as.matrix(tmp_cluster))
        }
      }
      fcs_include[[3]] <- as.matrix(add_cluster)
    }
  } else {
    fcs_include[[3]] <- NA
  }
  check_na <- sapply(fcs_include, function(arg1) {
    if(class(arg1)[1]=="logical") {
      if(is.na(arg1[1])) {
        return(NA)
      }
    } else {
      return(1)
    }
  })
  which_populated <- which(!is.na(check_na))
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
  if(subset_rows[1]!="all") {
    out_ff <- out_ff[subset_rows,]
  }
  if(include_timestamp) {
    flowCore::write.FCS(x = out_ff, filename = paste0(gsub("\\/$","",outdir),"/",gsub("\\.fcs$","",fcs_name),"_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".fcs"))
  } else {
    flowCore::write.FCS(x = out_ff, filename = paste0(gsub("\\/$","",outdir),"/",gsub("\\.fcs$","",fcs_name),".fcs"))
  }
}
