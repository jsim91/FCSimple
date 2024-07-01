fcs_batch_correction <- function(fcs_join_obj, correction_method = c("cyCombine"), correction_markers = "all",
                                 batch_source_regex = "[0-9]+\\-[A-Za-z]+\\-[0-9]+",
                                 cyCombine_SOMx = 8, cyCombine_SOMy = 8)
{
  # only cyCombine supported for now
  if(correction_method[1]=="cyCombine") {
    cyc_outdir1 <- paste0(getwd(),"/dbe_pre_correction")
    cyc_outdir2 <- paste0(getwd(),"/dbe_post_correction")
    rl <- readline(paste0("Based on getwd(), cyCombine detect_batch_effects will be saved here: ",cyc_outdir1,". Is this okay? (y/n) "))
    dir.create(path = cyc_outdir1, showWarnings = FALSE); dir.create(path = cyc_outdir2, showWarnings = FALSE)
    if(rl!="y") {
      stop("please update current working directory so that cyCombine::detect_batch_effects will be saved appropriately.")
    }
    require(cyCombine)
    require(stringr)
    exp_data <- fcs_join_obj$data
    if(correction_markers[1]=="all") {
      marks <- cyCombine::get_markers(df = exp_data)
    } else {
      marks <- correction_markers
    }
    exp_data$samples <- fcs_join_obj$source
    exp_data$batch <- stringr::str_extract(string = fcs_join_obj$source, pattern = batch_source_regex)
    print(paste0(length(unique(exp_data$batch))," batches found for correction: ",paste0(unique(exp_data$batch),collapse = ", ")))
    print("...detecting batch effects before correction...")
    cyCombine::detect_batch_effect(df = exp_data, out_dir = batch_correction_outdir, norm_method = "scale",
                                   xdim = cyCombine_SOMx, ydim = cyCombine_SOMy, seed = 123, batch_col = "batch",
                                   markers = marks)
    print(paste0("Consider reviewing cyCombine results here before proceeding beyond batch correction: ",cyc_outdir1))
    print("...normalizing data...")
    norm_data <- cyCombine::normalize(df = exp_data, markers = marks, norm_method = "scale")
    print("...creating som...")
    labels <- cyCombine::create_som(df = norm_data, markers = marks, rlen = 10, seed = 123)
    print("...correcting data... with anchor = NULL")
    corrected <- cyCombine::correct_data(df = exp_data, label = labels, markers = marks, covar = "batch")
    print("...deteching batch effects after correction...")
    cyCombine::detect_batch_effect(df = corrected, out_dir = cyc_outdir2, norm_method = "scale",
                                   xdim = cyCombine_SOMx, ydim = cyCombine_SOMy, seed = 123, batch_col = "batch",
                                   markers = marks)
    corrected_val <- as.data.frame(corrected)
    exprs_data_corrected <- corrected_val[,marks]
    corrected_source <- corrected_val$sample
    print("inserting corrected data into $data slot, uncorrected data will be stored under $data_uncorrected")
    print("saving correction history in $batch_correction_history")
    fcs_join_obj[['batch_correction_history']] <- list(method = "cyCombine", list(markers_corrected = marks,
                                                                                  batches_corrected = unique(exp_data$batch),
                                                                                  norm_method = "scale",
                                                                                  seed = 123,
                                                                                  SOMx_dim = cyCombine_SOMx,
                                                                                  SOMy_dim = cyCombine_SOMy,
                                                                                  datetime = Sys.time(),
                                                                                  session_info = sessionInfo()))
    fcs_join_obj[['data_uncorrected']] <- fcs_join_obj$data
    fcs_join_obj[['data']] <- exprs_data_corrected
    fcs_join_obj[['source']] <- corrected_source
    return(fcs_join_obj)
  }
}
