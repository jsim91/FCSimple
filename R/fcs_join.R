fcs_join <- function(files,
                     flowjo_diagnostics_file = NULL,
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
                     downsample_size = c(NA,25000),
                     batch_pattern = "[0-9]+\\-[A-Za-z]+\\-[0-9]+") {
  require(flowCore)
  if(any(length(files)==0,class(files[1])!="character")) {
    stop("'files' should be a vector of file names of .fcs files to be used in the analysis")
  }
  if(is.null(flowjo_diagnostics_file)) {
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
  }
  if(length(x = grep(pattern = "\\.fcs$", x = files, ignore.case = TRUE))!=length(files)) {
    files <- files[grep(pattern = "\\.fcs$", x = files, ignore.case = TRUE)]
    if(length(files)==0) {
      stop("error in argument 'files': No files with extension 'fcs' or 'FCS' found")
    }
  }
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
  if(use_descriptive_column_names) {
    colnames(raw_data) <- fs[[1]]@parameters@data$desc
  }
  if(!apply_transform) {
    return(list(data = raw_data,
                raw = raw_data,
                source = rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow))),
                object_history = paste0("joined: ",Sys.time())))
  }
  if(!is.null(flowjo_diagnostics_file)) {
    if(!grepl(pattern = "\\.txt$", x = flowjo_diagnostics_file)) {
      stop("'flowjo_diagnostics_file' must be the file name of a .txt file with contents copied from flowjo: Configure -> Diagnostics -> 'Show XML for Workspace'")
    } else {
      infcs <- fs[[1]]

      workspace <- read.delim(file = flowjo_diagnostics_file, header = FALSE, sep = "\n")
      tform_block_start <- grep(pattern = "<Transformations>", x = workspace[,1])[1]
      tform_block_end <- grep(pattern = "</Transformations>", x = workspace[,1])[1]
      transform_data <- workspace[tform_block_start:tform_block_end,1]
      desc_row <- grep(pattern = 'P[0-9]+(S|R|N)', workspace[,1])#[1:ncol(csv_list[[1]])]
      desc_row <- desc_row[-which(duplicated(workspace[desc_row,1]))]
      max_n <- max(as.numeric(gsub("P","",stringr::str_extract(string = workspace[desc_row,], pattern = "P[0-9]+"))))
      keyword_list <- vector("list", length = max_n); names(keyword_list) <- paste0("P",1:length(keyword_list))
      for(i in 1:length(keyword_list)) {
        tmpsubset <- workspace[desc_row,1][grep(pattern = paste0(names(keyword_list)[i],"(R|S|N)"), workspace[desc_row,1])]
        name_subset <- tmpsubset[grep("P[0-9]+N",tmpsubset)]; desc_subset <- tmpsubset[grep("P[0-9]+S.+value\\=[A-Za-z0-9]",tmpsubset)]
        name <- gsub("^ +|<Keyword name\\=| />|P[0-9]+N  value\\=|\\$","",name_subset)
        if(length(name)==0) {
          name <- ""
        }
        desc <- gsub("^ +|<Keyword name\\=| />|P[0-9]+S  value\\=|\\$","",desc_subset)
        if(length(desc)==0) {
          desc <- ""
        }
        if(name=="" & desc=="") {
          break
        }
        keyword_list[[i]] <- data.frame(desc = desc, name = name)
      }
      keyword_list <- keyword_list[which(sapply(keyword_list,class)=="data.frame")]
      param_data <- do.call(rbind,keyword_list); param_data$desc <- toupper(param_data$desc)

      block_starts <- grep(pattern = "<Transformations>", workspace[,1])
      block_ends <- grep(pattern = "</Transformations>", workspace[,1])
      transformations <- workspace[block_starts[1]:block_ends[1],1]

      tf1 <- grep(pattern = "<transforms:", x = transformations); tf2 <- grep(pattern = "</transforms:", x = transformations)
      if(length(tf1)!=length(tf2)) {
        stop("mismatch in index lengths..")
      }
      tf_list <- vector("list", length = length(tf1))
      if(use_fun=="hyperlog") {
        capture_dir <- system.file(package = "FCSimple")
      }
      for(i in 1:length(tf_list)) {
        tr_subset <- transformations[tf1[i]:tf2[i]]
        transform_fun_str <- stringr::str_extract(string = gsub("(^[ ]+<transforms:)", "", tr_subset[1]), pattern = "^[A-Za-z]+")
        spl_tr <- as.character(sapply(X = tr_subset[1], FUN = function(x) gsub("<|>|transforms:","",x)))
        spl_tr <- strsplit(x = spl_tr, split = " ")[[1]]; spl_tr <- spl_tr[which(spl_tr!="")]
        tf_list[[i]] <- list(fun = spl_tr[grep(pattern = "(linear|hyperlog|fasinh|biex)", spl_tr)], param = spl_tr[-grep(pattern = "(linear|hyperlog|fasinh|biex)", spl_tr)])
        elm2 <- strsplit(x = tr_subset[2], split = "=")[[1]]; elm2 <- gsub(" />","",elm2[length(elm2)])
        names(tf_list)[i] <- elm2
      }
      full_panel <- infcs@parameters@data
      tmp_raw <- raw_data
      tmp_data <- raw_data
      for(j in 1:ncol(tmp_data)) {
        if(!colnames(tmp_data)[j] %in% param_data$desc) {
          next
        }
        hyperparams <- tf_list[[param_data$name[which(param_data$desc==colnames(tmp_data)[j])[1]]]][[2]]
        use_fun <- tf_list[[param_data$name[which(param_data$desc==colnames(tmp_data)[j])[1]]]][[1]]
        if(!use_fun %in% c("hyperlog","biex","fasinh")) {
          stop(paste0("found function [",use_fun,"]: transformation using the flowjo workspace diagnostics file supports 'hyperlog', 'biex', and 'fasinh' transforms."))
        }
        if(use_fun=="hyperlog") {
          print(paste0("using hyperlog for ",colnames(tmp_data)[j]))
          hyper_t = as.numeric(gsub("T=","",hyperparams[grep("T=",hyperparams)]))
          hyper_a = as.numeric(gsub("A=","",hyperparams[grep("A=",hyperparams)]))
          hyper_m = as.numeric(gsub("M=","",hyperparams[grep("M=",hyperparams)]))
          hyper_w = as.numeric(gsub("W=","",hyperparams[grep("W=",hyperparams)]))
          if(hyper_a > (hyper_m - (2*hyper_w))) {
            hyper_a <- hyper_m - (2*hyper_w)
          }
          write.csv(x = tmp_data[,j], file = paste0(capture_dir,"/temp_files/__python_hyp_df__.csv"))
          system(command = paste0("python ",paste0(capture_dir,"/python/run_cluster.py")," ",
                                  paste0(capture_dir,"/temp_files/__python_cl_input__.mtx")," ",capture_dir,"/temp_files ",tolower(algorithm)," ",leiden_louvain_resolution))
          read_clus <- read.csv(paste0(capture_dir,"/temp_files/__tmp_cl__.csv"), check.names = FALSE)
          if(file.exists(paste0(capture_dir,"/temp_files/__tmp_cl__.csv"))) {
            file.remove(paste0(capture_dir,"/temp_files/__tmp_cl__.csv"))
          }
          if(file.exists(paste0(capture_dir,"/temp_files/__python_cl_input__.mtx"))) {
            file.remove(paste0(capture_dir,"/temp_files/__python_cl_input__.mtx"))
          }
          cluster_numbers <- read_clus[,1]
          # transform_fun <- flowCore::hyperlogtGml2(parameters = as.character(colnames(tmp_data)[j]),
          #                                          'T' = hyper_t,
          #                                          M = hyper_m,
          #                                          W = hyper_w,
          #                                          A = hyper_a)
          # tmp_data[,j] <- eval(transform_fun)(tmp_data)
          # some_data <- eval(transform_fun)(tmp_data)
          return(some_data)
        } else if(use_fun=="biex") {
          print(paste0("using biexp for ",colnames(tmp_data)[j]))
          widb <- as.numeric(gsub("width=","",hyperparams[grep("width=",hyperparams)]))
          transform_fun <- flowWorkspace::flowjo_biexp(pos = as.numeric(gsub("pos=","",hyperparams[grep("pos=",hyperparams)])),
                                                       neg = as.numeric(gsub("neg=","",hyperparams[grep("neg=",hyperparams)])),
                                                       # widthBasis = log(abs(widb),10)*-1)
                                                       widthBasis = as.numeric(gsub("width=","",hyperparams[grep("width=",hyperparams)])))
          tmp_data[,j] <- transform_fun(tmp_data[,j])
        } else if(use_fun=="fasinh") {
          print(paste0("using fasinh for ",colnames(tmp_data)[j]))
          transform_fun <- flowjo_fasinh(m = as.numeric(gsub("M=","",hyperparams[grep("M=",hyperparams)])),
                                         t = as.numeric(gsub("T=","",hyperparams[grep("T=",hyperparams)])),
                                         a = as.numeric(gsub("A=","",hyperparams[grep("A=",hyperparams)])))
          tmp_data[,j] <- transform_fun(tmp_data[,j])
        }
      }
      print("transformation completed successfully")
      src <- rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))
      return(list(data = tmp_data,
                  raw = tmp_raw,
                  source = src,
                  run_date = stringr::str_extract(string = src, pattern = batch_pattern),
                  transform_list = tf_list))
    }
  }
  if(!transform_per_channel) {
    if(tolower(instrument_type)=="cytof") {
      if(is.null(asinh_transform_cofactor[1])) {
        asinh_transform_cofactor <- 5
      }
      transform_function <- transformList(flowCore::colnames(fs), function(x) return(asinh(x/asinh_transform_cofactor)))
      fst <- flowCore::transform(fs, transform_function)
      for(i in 1:length(fst)) {
        if(i==1) {
          tmp_data <- flowCore::exprs(fst[[i]])
        } else {
          tmp_data <- rbind(tmp_data, flowCore::exprs(fst[[i]]))
        }
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
      src <- rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))
      return(list(data = tmp_data,
                  raw = raw_data,
                  source = src,
                  run_date = stringr::str_extract(string = src, pattern = batch_pattern),
                  object_history = paste0("joined: ",Sys.time())))
    } else if(tolower(instrument_type)=="flow") {
      if(transform_type=="asinh") {
        if(is.null(asinh_transform_cofactor[1])) {
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
      # if(length(grep("DATE|date|Date",names(fs[[1]]@description)))!=0) {
        src <- rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))
        return(list(data = tmp_data,
                    raw = raw_data,
                    source = src,
                    run_date = stringr::str_extract(string = src, pattern = batch_pattern),
                    object_history = paste0("joined: ",Sys.time())))
      # } else {
      #   return(list(data = tmp_data,
      #               raw = raw_data,
      #               source = rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow))),
      #               object_history = paste0("joined: ",Sys.time())))
      # }
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
                          source = rep(x = flowCore::sampleNames(fs),times = as.numeric(flowCore::fsApply(fs,nrow)))),
            file = paste0(system.file(package = "FCSimple"),"/temp_files/tmp_list_obj.rds"))
    require(shiny)
    shiny::runApp(appDir = file.path(system.file(package = "FCSimple"), "transform_app"))
  }
}
