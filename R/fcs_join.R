fcs_join <- function(files,
                     apply_transform = TRUE,
                     instrument_type = c("cytof","flow"),
                     use_descriptive_column_names = TRUE,
                     transform_function = NULL,
                     transform_type = c("asinh","biexp","hyperlog"),
                     asinh_transform_cofactor = 5,
                     biexp_transform_pos = 4.5,
                     biexp_transform_neg = 0,
                     biexp_transform_width = -10,
                     hyperlog_transform_T = 100000,
                     hyperlog_transform_M = 5,
                     hyperlog_transform_W = 0.001,
                     hyperlog_transform_A = 2,
                     transform_per_channel = TRUE,
                     downsample_size = c(NA,25000)) {
  require(flowCore)
  if(!transform_per_channel) {
    if(length(instrument_type)>1) {
      warning(paste0("Consider specifying 'instrument_type'. Default use is 'cytof'. If inputs are from a flow cytometer, use 'flow'. Using ",instrument_type[1]," for now."))
      instrument_type <- instrument_type[1]
    }
    if(length(transform_type)>1) {
      warning(paste0("Consider specifying 'transform_type'. Default is 'asinh'. If undesirable, you can specify 'biexp' or 'hyperlog'. Using ",transform_type[1]," for now."))
      transform_type <- transform_type[1]
    }
  }
  if(length(x = grep(pattern = "\\.fcs$", x = files, ignore.case = TRUE))!=length(files)) {
    files <- files[grep(pattern = "\\.fcs$", x = files, ignore.case = TRUE)]
    if(length(files)==0) {
      stop("error in argument 'files': No files with extension 'fcs' or 'FCS' found")
    }
  }
  # fcs_list <- vector("list", length = length(files)); names(fcs_list) <- files
  # for(i in 1:length(fcs_list)) {
  #   tmp_fcs <- read.FCS(filename = files[i], truncate_max_range = FALSE)
  #   na_channel <- which(is.na(tmp_fcs@parameters@data$desc))
  #   if(length(na_channel)!=0) {
  #     tmp_fcs <- tmp_fcs[,-na_channel]
  #   }
  #   if(!is.na(downsample_size)) {
  #     if(nrow(tmp_fcs)>downsample_size) {
  #       set.seed(123)
  #       fcs_list[[i]] <- tmp_fcs[sample(1:nrow(tmp_fcs),downsample_size,replace=FALSE),]
  #       colnames(fcs_list[[i]]) <- fcs_list[[1]]@parameters@data$desc
  #     } else {
  #       fcs_list[[i]] <- tmp_fcs
  #       colnames(fcs_list[[i]]) <- fcs_list[[1]]@parameters@data$desc
  #     }
  #   } else {
  #     fcs_list[[i]] <- tmp_fcs
  #     colnames(fcs_list[[i]]) <- fcs_list[[1]]@parameters@data$desc
  #   }
  # }
  # fs <- flowSet(fcs_list)
  fs <- flowCore::read.flowSet(files = files, truncate_max_range = FALSE)
  if(mean(is.na(downsample_size))!=0) {
    downsample_size <- NA
  } else if(length(downsample_size)!=0){
    downsample_size <- downsample_size[1]
  }
  sampleNames(fs) <- gsub("^.+/","",sampleNames(fs))
  if(!is.na(downsample_size)) {
    for(i in 1:length(fs)) {
      if(nrow(fs[[i]])>downsample_size) {
        set.seed(123)
        fs[[i]] <- fs[[i]][sample(x = 1:nrow(fs[[i]]), size = downsample_size, replace = FALSE),]
      }
    }
  }
  for(i in 1:length(fs)) {
    if(i==1) {
      raw_data <- exprs(fs[[i]])
    } else {
      raw_data <- rbind(raw_data, exprs(fs[[i]]))
    }
  }
  # if(use_ncdf) {
  #   require(ncdfFlow)
  # }
  # if(use_ncdf) {
  #   fs <- ncdfFlow::read.ncdfFlowSet(files = files)
  # } else {
  #   fs <- flowCore::read.flowSet(files = files, truncate_max_range = FALSE)
  # }
  if(!apply_transform) {
    return(list(data = raw_data,
                raw = raw_data,
                source = rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))))
  }
  if(!transform_per_channel) {
    if(tolower(instrument_type)=="cytof") {
      if(is.null(asinh_transform_cofactor)) {
        asinh_transform_cofactor <- 5
      }
      transform_function <- transformList(flowCore::colnames(fs), function(x) return(asinh(x/asinh_transform_cofactor)))
      # if(use_ncdf) {
      #   fst <- ncdfFlow::transform(fs, transform_function)
      # } else {
      fst <- flowCore::transform(fs, transform_function)
      # }
      for(i in 1:length(fst)) {
        if(i==1) {
          tmp_data <- flowCore::exprs(fst[[i]])
        } else {
          tmp_data <- rbind(tmp_data, flowCore::exprs(fst[[i]]))
        }
      }
      return(list(data = tmp_data,
                  raw = raw_data,
                  source = rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))))
    } else if(tolower(instrument_type)=="flow") {
      if(transform_type=="asinh") {
        if(is.null(asinh_transform_cofactor)) {
          asinh_transform_cofactor <- 200
        }
        transform_function <- flowCore::transformList(flowCore::colnames(fs), function(x) return(asinh(x/asinh_transform_cofactor)))
        fst <- flowCore::transform(fs, transform_function)
        for(i in 1:length(fst)) {
          if(i==1) {
            tmp_data <- flowCore::exprs(fst[[i]])
          } else {
            tmp_data <- rbind(tmp_data, flowCore::exprs(fst[[i]]))
          }
        }
      } else if(transform_type=="biexp") {
        require(flowWorkspace)
        if(any(!is.numeric(biexp_transform_pos), !is.numeric(biexp_transform_neg), !is.numeric(biexp_transform_width))) {
          stop("error in argument(s) 'biexp_transform_.': values must be numeric")
        }
        transform_function <- flowWorkspace::flowjo_biexp(pos = biexp_transform_pos,
                                                          neg = biexp_transform_neg,
                                                          widthBasis = biexp_transform_width)
        for(i in 1:length(fs)) {
          if(i==1) {
            tmp_data <- flowCore::exprs(object = fs[[i]])
          } else {
            tmp_data <- rbind(tmp_data, flowCore::exprs(object = fs[[i]]))
          }
        }
        for(i in 1:ncol(tmp_data)) {
          tmp_data[,i] <- transform_function(tmp_data[,i])
        }
      } else if(transform_type=="hyperlog") {
        stop("Run this function with 'apply_transform = FALSE' then use fcs_as.hyperlog to transform the values (in $data slot) using the hyperlog transform.")
        # if(any(!is.numeric(hyperlog_transform_T), !is.numeric(hyperlog_transform_M), !is.numeric(hyperlog_transform_W), !is.numeric(hyperlog_transform_A))) {
        #   stop("error in argument(s) 'hyperlog_transform_.': values must be numeric")
        # }
        # transf_hyperlog <- function(fset, hyper_t, hyper_m, hyper_w, hyper_a) {
        #   for(i in 1:length(fset)) {
        #     exprs_data <- exprs(fset[[i]])
        #     for(j in 1:ncol(exprs_data)) {
        #       transform_fun <- flowCore::hyperlogtGml2(parameters = colnames(exprs_data)[j], 'T' = hyper_t, M = hyper_m, W = hyper_w, A = hyper_a, transformationId = "hyper1")
        #       exprs_data[,j] <- eval(transform_fun)(exprs_data)
        #     }
        #     if(i==1) {
        #       transf_data <- exprs_data
        #     } else {
        #       transf_data <- rbind(transf_data, exprs_data)
        #     }
        #   }
        #   return(transf_data)
        # }
        # tmp_data <- transf_hyperlog(fset = fs, hyper_t = hyperlog_transform_T, hyper_m = hyperlog_transform_M,
        #                             hyper_w = hyperlog_transform_W, hyper_a = hyperlog_transform_A)
      } else {
        stop("error in argument 'instrument_type': depending on how FCS files were created, use 'cytof' or 'flow'")
      }
      if(nrow(tmp_data)==0) {
        stop("error: data set contains 0 events (rows)")
      }
      if(ncol(tmp_data)==0) {
        stop("error: data set contains 0 channels (columns)")
      }
      if(use_descriptive_column_names) {
        desc_names <- fs[[1]]@parameters@data$desc
        if(length(desc_names)==ncol(tmp_data)) {
          colnames(tmp_data) <- desc_names
          colnames(raw_data) <- desc_names
        } else {
          print("Unable to find descriptive column names. Using original names.")
        }
      }
      if(length(grep("DATE|date|Date",names(fs[[1]]@description)))!=0) {
        run_dates <- flowCore::fsApply(fs, function(x) return(x@description[[grep("DATE|date|Date",names(x@description))[1]]]))
        return(list(data = tmp_data,
                    raw = raw_data,
                    source = rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow))),
                    run_date = ifelse(length(run_dates)>0,run_dates,NULL)))
      } else {
        return(list(data = tmp_data,
                    raw = raw_data,
                    source = rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))))
      }
    }
  } else {
    for(i in 1:length(fs)) {
      if(i==1) {
        tmp_data <- flowCore::exprs(object = fs[[i]])
      } else {
        tmp_data <- rbind(tmp_data, flowCore::exprs(object = fs[[i]]))
      }
    }
    if(use_descriptive_column_names) {
      desc_names <- fs[[1]]@parameters@data$desc
      if(length(desc_names)==ncol(tmp_data)) {
        colnames(tmp_data) <- desc_names
        colnames(raw_data) <- desc_names
      } else {
        colnames(tmp_data) <- fs[[1]]@parameters@data$name
        print("Unable to find descriptive column names. Using original names.")
      }
    }
    temp_files <- list.files(path = paste0(system.file(package = "FCSimple"),"/temp_files/"), full.names = TRUE, recursive = TRUE)
    if(length(temp_files)!=0) {
      file.remove(temp_files)
    }
    if(nrow(tmp_data)>50000) {
      set.seed(123)
      write.csv(x = tmp_data[sample(1:nrow(tmp_data),size=50000,replace=F),],
                file = paste0(system.file(package = "FCSimple"),"/temp_files/tmp_data.csv"), row.names = FALSE)
    } else {
      write.csv(x = tmp_data,file = paste0(system.file(package = "FCSimple"),"/temp_files/tmp_data.csv"), row.names = FALSE)
    }
    saveRDS(object = list(data = tmp_data,
                          source = rep(x = flowCore::sampleNames(fs),
                                       times = as.numeric(flowCore::fsApply(fs,nrow)))),
            file = paste0(system.file(package = "FCSimple"),"/temp_files/tmp_list_obj.rds"))
    require(shiny)
    shiny::runApp(appDir = file.path(system.file(package = "FCSimple"), "transform_app"))
  }
}
