fcs_write.FCS <- function(fcs_join_obj,
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
    for(i in 1:length(include_reductions)) {
      if(tolower(include_reductions[i]) %in% names(fcs_join_obj)) {
        append_coords <- as.matrix(fcs_join_obj[[tolower(include_reductions[i])]][["coordinates"]])
        if(i==1) {
          data_return_1 <- append_coords
        } else {
          data_return_1 <- cbind(data_return_1, append_coords)
        }
      } else {
        warning(paste0("results of reduction ",tolower(include_reductions[i])," not found in object and will not be included"))
      }
    }
    fcs_include[[2]] <- data_return_1
  }
  if(length(include_clusterings)!=0) {
    cluster_numbers <- vector("list", length = length(include_clusterings))
    names(cluster_numbers) <- paste0(tolower(include_clusterings),"_cluster")
    for(i in 1:length(include_clusterings)) {
      if(tolower(include_clusterings[i]) %in% names(fcs_join_obj)) {
        cluster_numbers[[i]] <- fcs_join_obj[[tolower(include_clusterings[i])]][["clusters"]]
      } else {
        warning(paste0("results of clustering ",tolower(include_clusterings[i])," not found in object and will not be included"))
      }
    }
    rm_elements <- c()
    for(i in 1:length(cluster_numbers)) {
      if(length(cluster_numbers[[i]])==0) {
        rm_elements <- append(rm_elements,i)
      }
    }
    if(length(rm_elements)!=0) {
      if(length(rm_elements)==1) {
        cluster_numbers[[rm_elements]] <- NULL
      } else {
        cluster_numbers[rm_elements] <- NULL
      }
    }
    for(i in 1:length(cluster_numbers)) {
      append_numbers <- as.matrix(data.frame(var1 = cluster_numbers[[i]]))
      colnames(append_numbers) <- names(cluster_numbers)[i]
      if(i==1) {
        data_return_2 <- append_numbers
      } else {
        data_return_2 <- cbind(data_return_2, append_numbers)
      }
    }
    fcs_include[[3]] <- data_return_2
  }
  which_populated <- which(sapply(fcs_include, function(x) return(class(x)[1]))!="NULL")
  if(length(which_populated)<3) {
    if(length(which_populated)==1) {
      final_return <- fcs_include[[which_populated]]
    } else if(length(which_populated)==2) {
      fcs_include[which(!1:3 %in% which_populated)] <- NULL
      final_return <- cbind(fcs_include[[1]], fcs_include[[2]])
    }
  } else if(length(which_populated)==3) {
    for(i in 1:length(which_populated)) {
      if(i==1) {
        next
      } else if(i==2) {
        final_return <- cbind(fcs_include[[1]], fcs_include[[i]])
      } else {
        final_return <- cbind(final_return, fcs_include[[i]])
      }
    }
  }
  if(class(final_return)[1]!="matrix") {
    final_return <- as.matrix(final_return)
  }
  out_ff <- new("flowFrame", exprs = final_return)
  if(include_timestamp) {
    write.FCS(x = out_ff, filename = paste0(gsub("\\/$","",outdir,"/fcsimple_outfile_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".fcs")))
  } else {
    write.FCS(x = out_ff, filename = paste0(gsub("\\/$","",outdir,"/fcsimple_outfile.fcs")))
  }
}
