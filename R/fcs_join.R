#' @title Read and Join Multiple FCS or CSV Files into a Single Analysis Object
#'
#' @description
#'   Reads a set of FCS or CSV files, optionally downsamples and applies transformations
#'   (FlowJo diagnostics, asinh, biexp, or hyperlog), and concatenates them into
#'   a unified data matrix with accompanying metadata and history.
#'   
#'   For CSV files, the function assumes data is already transformed. CSV files should 
#'   be exported from FlowJo using 'CSV - Channel values' with 'Include header' and 
#'   'Stain' selected after all transforms have been applied. Ensure a run date ($DATE) 
#'   is included in the file names.
#'
#' @param files
#'   Character vector of file paths to .fcs or .csv files. For CSV files, data 
#'   should be pre-transformed and exported from FlowJo with proper formatting 
#'   (see Details).
#'
#' @param flowjo_diagnostics_file
#'   Optional path to a FlowJo diagnostics TXT file (Configure → Diagnostics →
#'   “Show XML for Workspace”). Provides per-channel transform parameters.
#'   Default `NULL`.
#'
#' @param apply_transform
#'   Logical; if `TRUE` (default), applies channel transformations. If `FALSE`,
#'   returns raw data only.
#'
#' @param instrument_type
#'   Deprecated. Instrument type is now auto-detected: if min(data) >= 0, assumes
#'   'cytof' (mass cytometry); otherwise assumes 'flow' (fluorescence cytometry).
#'
#' @param use_descriptive_column_names
#'   Logical; if `TRUE` (default), replaces channel names with descriptive labels
#'   from FCS metadata.
#'
#' @param transform_function
#'   Optional user-supplied `transformList` (e.g., from
#'   `flowCore::transformList`) for per-channel transforms. Overrides
#'   `transform_type` when `transform_per_channel = TRUE`.
#'
#' @param transform_type
#'   Character; global transform to apply when `transform_per_channel = FALSE`.
#'   One of `"asinh"`, `"biexp"`, or `"hyperlog"`. Default `"asinh"`.
#'
#' @param asinh_transform_cofactor
#'   Numeric cofactor for asinh transforms (default 5).
#'
#' @param biexp_transform_pos
#'   Numeric positive breakpoint for biexp (default 4.5).
#'
#' @param biexp_transform_neg
#'   Numeric negative breakpoint for biexp (default 0).
#'
#' @param biexp_transform_width
#'   Numeric width basis for biexp (default -10).
#'
#' @param hyperlog_transform_T
#'   Numeric “T” parameter for hyperlog transform (default 100000).
#'
#' @param hyperlog_transform_M
#'   Numeric “M” decades for hyperlog transform (default 5).
#'
#' @param hyperlog_transform_W
#'   Numeric “W” width around zero for hyperlog transform (default 0.001).
#'
#' @param hyperlog_transform_A
#'   Numeric “A” offset for hyperlog transform (default 2).
#'
#' @param transform_per_channel
#'   Logical; if `TRUE` (default), launches an interactive Shiny application
#'   (in `transform_app/`) to preview and tweak per-channel transforms. The
#'   function will block until you finish and close the app. If `FALSE`,
#'   applies the global `transform_type`.
#'
#' @param downsample_size
#'   Integer or `NA`; maximum events per file to sample before joining.
#'   Default `NA` (no downsampling).
#'
#' @param batch_pattern
#'   Regular expression to extract `run_date` from sample names in `source`.
#'   Default `"[0-9]+\\-[A-Za-z]+\\-[0-9]+"`.
#'
#' @details
#'   Internally, the function:
#'   1. For FCS files: Reads all files into a `flowSet` via `flowCore::read.flowSet()`.
#'      For CSV files: Reads using `data.table::fread()` assuming pre-transformed data.
#'   2. Optionally downsamples each file to `downsample_size`.
#'   3. Extracts expression matrices with `flowCore::exprs()` (FCS) or directly from 
#'      data frames (CSV).
#'   4. If `transform_per_channel = TRUE` (FCS only), launches the package's
#'      `transform_app/` Shiny GUI so you can inspect per-channel
#'      diagnostics and build a custom `transformList`.
#'   5. Otherwise, applies the global transform as specified by
#'      `transform_type` or user's `transform_function` (FCS only).
#'   6. Concatenates raw and transformed data into `raw` and `data` matrices.
#'   7. Builds `source` and `run_date` via `flowCore::sampleNames()` or file names 
#'      and `stringr::str_extract()`.
#'   8. Creates `metadata` data frame with patient_ID and run_date.
#'   9. Records `collection_instrument` and appends a timestamp in
#'      `object_history`.
#'
#' @return
#'   A list with elements:
#'   - `data`: numeric matrix of transformed events × channels (or data frame for CSV).
#'   - `raw`: numeric matrix of untransformed events × channels (or NA for CSV files).
#'   - `source`: character vector of sample identifiers for each event.
#'   - `run_date`: character vector parsed from `source` via `batch_pattern`.
#'   - `metadata`: data frame with columns `patient_ID` and `run_date` (one row per file).
#'   - `transform_list`: per-channel transform parameters (if FlowJo diagnostics file used).
#'   - `collection_instrument`: chosen `instrument_type`.
#'   - `object_history`: timestamped entry recording the join.
#'
#' @examples
#' \dontrun{ (FCS files)
#' files <- list.files("data/fcs", pattern = "\\.fcs$", full.names = TRUE)
#' joined <- FCSimple::fcs_join(files)
#'
#' # Interactive per-channel diagnostics via Shiny
#' joined2 <- FCSimple::fcs_join(
#'   files,
#'   transform_per_channel = TRUE
#' )
#'
#' # Raw join without any transforms
#' joined_raw <- FCSimple::fcs_join(files, apply_transform = FALSE)
#'
#' # Join pre-transformed CSV files exported from FlowJo
#' csv_files <- list.files("data/flowjo_export", pattern = "\\.csv$", full.names = TRUE)
#' joined_csv <- FCSimple::fcs_join(csv_files
#' joined_raw <- FCSimple::fcs_join(files, apply_transform = FALSE)
#' }
#'
#' @seealso
#'   flowCore::read.flowSet, flowCore::exprs, flowCore::transformList,
#'   flowWorkspace::flowjo_biexp, stringr::str_extract
#'
#' @importFrom flowCore read.flowSet exprs sampleNames fsApply transformList
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom stringr str_extract
#' @export
fcs_join <- function(files,
                     flowjo_diagnostics_file = NULL,
                     apply_transform = TRUE,
                     use_descriptive_column_names = TRUE,
                     transform_function = NULL,
                     transform_type = "asinh",
                     asinh_transform_cofactor = 5,
                     biexp_transform_pos = 4.5,
                     biexp_transform_neg = 0,
                     biexp_transform_width = -10,
                     hyperlog_transform_T = 100000,
                     hyperlog_transform_M = 5,
                     hyperlog_transform_W = 0.001,
                     hyperlog_transform_A = 2,
                     transform_per_channel = TRUE,
                     downsample_size = NA,
                     batch_pattern = "[0-9]+\\-[A-Za-z]+\\-[0-9]+") {
  # testing
  # files <- list.files(path = "J:/oakes_flow/scenith_full_pilot/box_dl", full.names = TRUE)
  # flowjo_diagnostics_file = "J:/oakes_flow/scenith_full_pilot/josh_transforms/workspace_transforms_2.txt"
  # apply_transform = TRUE
  # instrument_type = "flow"
  # use_descriptive_column_names = TRUE
  # transform_function = NULL
  # transform_type = "hyperlog"
  # asinh_transform_cofactor = 5
  # biexp_transform_pos = 4.5
  # biexp_transform_neg = 0
  # biexp_transform_width = -10
  # hyperlog_transform_T = 100000
  # hyperlog_transform_M = 5
  # hyperlog_transform_W = 0.001
  # hyperlog_transform_A = 2
  # transform_per_channel = TRUE
  # downsample_size = NA
  # batch_pattern = "[0-9]+\\-[A-Za-z]+\\-[0-9]+"


  require(flowCore)
  # oo <- options(scipen = 100000000000)
  # on.exit(options(oo))
  options(scipen = 1000000)
  if(any(length(files)==0,class(files[1])!="character")) {
    stop("'files' should be a vector of file names of .fcs or .csv files to be used in the analysis")
  }
  if(mean(grepl(pattern = '\\.csv$', x = files))==1) { # csv case
    warning("'fcs_join with csv files assumes data is already transformed. Export from flowjo using 'CSV - Channel values' with 'Include header' and 'Stain' selected after all transforms have been sent. Ensure a run date ($DATE) is included in the file names. Refer to https://github.com/jsim91/FCSimple 'Data Input Format' section.")
    csv_data <- lapply(X = files, FUN = data.table::fread, check.names = FALSE,
                       data.table = FALSE, nThread = parallel::detectCores())
    names(csv_data) <- gsub(pattern = '^.\\/', replacement = '', x = files)
    if(!is.na(downsample_size)) {
      csv_ds <- lapply(X = csv_data, FUN = function(arg) {
        if(nrow(arg)>downsample_size) {
          set.seed(123); return(arg[sample(1:nrow(arg), downsample_size, replace = F),])
        } else {
          return(arg)
        }
      })
    } else {
      csv_ds <- csv_data
    }
    rm(csv_data); gc()
    csv <- do.call(rbind, csv_data)
    src <- rep(names(csv_data), rep(sapply(csv_data, nrow)))
    rd <- stringr::str_extract(src, batch_pattern)
    meta <- data.frame(patient_ID = src, run_date = rd); meta <- meta[!duplicated(meta$patient_ID),]
    
    # Auto-detect instrument type based on data
    detected_instrument <- if(min(csv, na.rm = TRUE) >= 0) "cytof" else "flow"
    
    return(list(data = csv_ds,
                raw = NA,
                source = src,
                run_date = rd,
                metadata = meta,
                collection_instrument = detected_instrument,
                object_history = paste0("joined: ",Sys.time())))
  }
  if(is.null(flowjo_diagnostics_file)) {
    if(!transform_per_channel) {
      # No warnings needed - defaults are now sensible
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
      raw_data <- flowCore::exprs(fs[[i]])
    } else {
      raw_data <- rbind(raw_data, flowCore::exprs(fs[[i]]))
    }
  }
  if(use_descriptive_column_names) {
    colnames(raw_data) <- fs[[1]]@parameters@data$desc
  }
  if(!apply_transform) {
    src <- rep(x = flowCore::sampleNames(fs), times = as.numeric(flowCore::fsApply(fs,nrow)))
    rd <- stringr::str_extract(string = src, pattern = batch_pattern)
    if(sum(is.na(rd))!=0) {
      warning('One or more run_date is NA. Returning placeholder run_date values. Consider adjusting batch_pattern.')
      rd <- rep('placeholder')
    }
    meta <- data.frame(patient_ID = src, run_date = rd); meta <- meta[!duplicated(meta$patient_ID),]
    
    # Auto-detect instrument type based on data
    detected_instrument <- if(min(raw_data, na.rm = TRUE) >= 0) "cytof" else "flow"
    
    return(list(data = raw_data,
                raw = raw_data,
                source = src,
                run_date = rd,
                metadata = meta,
                collection_instrument = detected_instrument,
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
      which_dupl <- which(duplicated(workspace[desc_row,1]))
      if(length(which_dupl)!=0) {
        desc_row <- desc_row[-which(duplicated(workspace[desc_row,1]))]
      }
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
        if(name=="" && desc=="") {
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
        if(!toupper(colnames(tmp_data)[j]) %in% toupper(param_data$desc)) {
          next
        }
        hyperparams <- tf_list[[param_data$name[which(toupper(param_data$desc)==toupper(colnames(tmp_data)[j]))[1]]]][[2]]
        use_fun <- tf_list[[param_data$name[which(toupper(param_data$desc)==toupper(colnames(tmp_data)[j]))[1]]]][[1]]
        if(!use_fun %in% c("hyperlog","biex","fasinh")) {
          stop(paste0("found function [",use_fun,"]: transformation using the flowjo workspace diagnostics file supports 'hyperlog', 'biex', and 'fasinh' transforms."))
        }
        if(use_fun=="hyperlog") {
          print(paste0("using hyperlog for ",colnames(tmp_data)[j]))
          capture_dir <- system.file(package = "FCSimple")
          hyper_t = as.numeric(gsub("T=","",hyperparams[grep("T=",hyperparams)]))
          hyper_w = as.numeric(gsub("W=","",hyperparams[grep("W=",hyperparams)]))
          hyper_m = as.numeric(gsub("M=","",hyperparams[grep("M=",hyperparams)]))
          hyper_a = as.numeric(gsub("A=","",hyperparams[grep("A=",hyperparams)]))
          if(hyper_a > (hyper_m - (2*hyper_w))) {
            hyper_a <- hyper_m - (2*hyper_w)
          }
          write.csv(x = tmp_data[,j], file = paste0(capture_dir,"/temp_files/__python_hyp_df__.csv"), row.names = FALSE)
          system(command = paste0("python ",paste0(capture_dir,"/python/transf_hyperlog.py")," ",paste0(capture_dir,"/temp_files/__python_hyp_df__.csv")," ",capture_dir,"/temp_files ",
                                  hyper_t," ",hyper_w," ",hyper_m," ",hyper_a))
          read_exprs <- read.csv(paste0(capture_dir,"/temp_files/__tmp_exprs__.csv"), check.names = FALSE)
          if(file.exists(paste0(capture_dir,"/temp_files/__tmp_exprs__.csv"))) {
            file.remove(paste0(capture_dir,"/temp_files/__tmp_exprs__.csv"))
          }
          if(file.exists(paste0(capture_dir,"/temp_files/__python_hyp_df__.csv"))) {
            file.remove(paste0(capture_dir,"/temp_files/__python_hyp_df__.csv"))
          }
          tmp_data[,j] <- read_exprs[,1]
          # transform_fun <- flowCore::hyperlogtGml2(parameters = as.character(colnames(tmp_data)[j]),
          #                                          'T' = hyper_t,
          #                                          M = hyper_m,
          #                                          W = hyper_w,
          #                                          A = hyper_a)
          # tmp_data[,j] <- eval(transform_fun)(tmp_data)
          # some_data <- eval(transform_fun)(tmp_data)
          # return(some_data)
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
      rd <- stringr::str_extract(string = src, pattern = batch_pattern)
      if(sum(is.na(rd))!=0) {
        warning('One or more run_date is NA. Returning placeholder run_date values. Consider adjusting batch_pattern.')
        rd <- rep('placeholder')
      }
      meta <- data.frame(patient_ID = src, run_date = rd); meta <- meta[!duplicated(meta$patient_ID),]
      
      # Auto-detect instrument type based on raw data
      detected_instrument <- if(min(tmp_raw, na.rm = TRUE) >= 0) "cytof" else "flow"
      
      return(list(data = tmp_data,
                  raw = tmp_raw,
                  source = src,
                  run_date = rd,
                  metadata = meta,
                  transform_list = tf_list,
                  collection_instrument = detected_instrument,
                  object_history = paste0("joined: ",Sys.time())))
    }
  }
  if(!transform_per_channel) {
    # Auto-detect instrument type from raw data to choose appropriate transform cofactor
    detected_instrument <- if(min(raw_data, na.rm = TRUE) >= 0) "cytof" else "flow"
    
    if(detected_instrument == "cytof") {
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
      rd <- stringr::str_extract(string = src, pattern = batch_pattern)
      if(sum(is.na(rd))!=0) {
        warning('One or more run_date is NA. Returning placeholder run_date values. Consider adjusting batch_pattern.')
        rd <- rep('placeholder')
      }
      meta <- data.frame(patient_ID = src, run_date = rd); meta <- meta[!duplicated(meta$patient_ID),]
      
      return(list(data = tmp_data,
                  raw = raw_data,
                  source = src,
                  run_date = rd,
                  metadata = meta,
                  collection_instrument = detected_instrument,
                  object_history = paste0("joined: ",Sys.time())))
    } else if(detected_instrument == "flow") {
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
        rd <- stringr::str_extract(string = src, pattern = batch_pattern)
        if(sum(is.na(rd))!=0) {
          warning('One or more run_date is NA. Returning placeholder run_date values. Consider adjusting batch_pattern.')
          rd <- rep('placeholder')
        }
        meta <- data.frame(patient_ID = src, run_date = rd); meta <- meta[!duplicated(meta$patient_ID),]
        
        # Auto-detect instrument type based on data
        detected_instrument <- if(min(tmp_data, na.rm = TRUE) >= 0) "cytof" else "flow"
        
        return(list(data = tmp_data,
                    raw = raw_data,
                    source = src,
                    run_date = rd,
                    metadata = meta,
                    collection_instrument = detected_instrument,
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
