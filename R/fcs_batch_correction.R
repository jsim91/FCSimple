fcs_batch_correction <- function(fcs_join_obj, use_rep = "data", correction_method = c("cyCombine", "harmony"),
                                 correction_markers = "all", batch_source_regex = "[0-9]+\\-[A-Za-z]+\\-[0-9]+",
                                 cyCombine_SOMx = 8, cyCombine_SOMy = 8, cyCombine_detect_effects = FALSE,
                                 harmony_cores = 1, harmony_iterations = 10, harmony_covars = c("batch","sample"),
                                 harmony_lambda = 1, harmony_sample_element = 'source')
{
  # larger harmony_lambda (ridge regression penalty) protect against over correction
  use_rep <- tolower(use_rep)
  if(!use_rep %in% c("data","pca")) {
    stop("'use_rep' indicates what representation of the data will be passed to the algorithm. Use 'data' or 'pca'. 'pca' requires that FCSimple::fcs_pca was already run on the data.")
  } else {
    if(use_rep=="data") {
      rep_data <- fcs_join_obj[["data"]]
      print("Using fcs_join_obj[['data']] for correction")
    } else if(use_rep=="pca") {
      rep_data <- fcs_join_obj[["pca"]][['pca_data']]
      print("Using fcs_join_obj[['pca']][['pca_data']] for correction")
    }
  }
  if(correction_method[1]=="cyCombine") {
    cmeth <- "cyCombine"
    require(cyCombine)
    require(stringr)
    exp_data <- as.data.frame(rep_data)
    if(correction_markers[1]=="all") {
      marks <- cyCombine::get_markers(df = exp_data)
    } else {
      marks <- correction_markers
    }
    # exp_data$samples <- fcs_join_obj$source
    exp_data$samples <- as.numeric(factor(fcs_join_obj$source))
    # exp_data$batch <- stringr::str_extract(string = fcs_join_obj$source, pattern = batch_source_regex)
    exp_data$batch <- as.numeric(factor(fcs_join_obj[["run_date"]]))
    print(paste0(length(unique(exp_data$batch))," batches found for correction: ",paste0(unique(exp_data$batch),collapse = ", ")))
    if(cyCombine_detect_effects) {
      cyc_outdir1 <- paste0(getwd(),"/dbe_pre_correction")
      cyc_outdir2 <- paste0(getwd(),"/dbe_post_correction")
      # rl <- readline(paste0("Based on getwd(), cyCombine detect_batch_effects will be saved here: ",cyc_outdir1,". Is this okay? (y/n) "))
      paste0("Based on getwd(), cyCombine detect_batch_effects will be saved here: ",cyc_outdir1," ... save location uses getwd() as root.")
      dir.create(path = cyc_outdir1, showWarnings = FALSE); dir.create(path = cyc_outdir2, showWarnings = FALSE)
      # if(rl!="y") {
      #   stop("please update current working directory so that cyCombine::detect_batch_effects will be saved appropriately.")
      # }
      require(outliers)
      print("...detecting batch effects before correction...")
      cyCombine::detect_batch_effect(df = exp_data, out_dir = cyc_outdir1, norm_method = "scale",
                                     xdim = cyCombine_SOMx, ydim = cyCombine_SOMy, seed = 123, batch_col = "batch",
                                     markers = marks)
      print(paste0("pre-batch correction: ",cyc_outdir1))
    }
    print("...normalizing data...")
    norm_data <- cyCombine::normalize(df = exp_data, markers = marks, norm_method = "scale")
    print("...creating som...")
    labels <- cyCombine::create_som(df = norm_data, markers = marks, rlen = 10, seed = 123)
    print("...correcting data... with anchor = NULL")
    corrected <- cyCombine::correct_data(df = exp_data, label = labels, markers = marks, covar = "batch")
    if(cyCombine_detect_effects) {
      print("...detecting batch effects after correction...")
      cyCombine::detect_batch_effect(df = corrected, out_dir = cyc_outdir2, norm_method = "scale",
                                     xdim = cyCombine_SOMx, ydim = cyCombine_SOMy, seed = 123, batch_col = "batch",
                                     markers = marks)
      print(paste0("post-batch correction: ",cyc_outdir2))
    }
    corrected_val <- as.data.frame(corrected)
    exprs_data_corrected <- corrected_val[,marks]
    corrected_source <- corrected_val$sample
    print("saving correction history in $batch_correction")
    fcs_join_obj[['batch_correction']] <- list(data = exprs_data_corrected,
                                               source = corrected_source,
                                               method = cmeth,
                                               other = list(markers_corrected = marks,
                                                            batches_corrected = unique(exp_data$batch),
                                                            norm_method = "scale",
                                                            seed = 123,
                                                            SOMx_dim = cyCombine_SOMx,
                                                            SOMy_dim = cyCombine_SOMy,
                                                            datetime = Sys.time(),
                                                            session_info = sessionInfo()))
  } else if(correction_method[1]=="harmony") {
    require(harmony)
    cmeth <- "harmony"
    harm_in <- as.matrix(rep_data)
    print(paste0("batches found for harmony correction: ", paste0(unique(fcs_join_obj[["run_date"]]), collapse = ", ")))
    harm_meta <- data.frame(cell_id = 1:nrow(harm_in), batch = fcs_join_obj[["run_date"]], sample = fcs_join_obj[[harmony_sample_element]])
    harm_out <- harmony::RunHarmony(data_mat = harm_in, meta_data = harm_meta, vars_use = harmony_covars, ncores = harmony_cores,
                                    max_iter = harmony_iterations, lambda = harmony_lambda)
    fcs_join_obj[['batch_correction']] <- list(data = harm_out,
                                               harmony_meta = harm_meta,
                                               method = cmeth,
                                               other = list(number_of_batches = length(unique(harm_meta$batch)),
                                                            batches = unique(harm_meta$batch),
                                                            datetime = Sys.time(),
                                                            session_info = sessionInfo()))
  }
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  } else {
    fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(cmeth," on ",use_rep,": ",Sys.time()))
  }
  return(fcs_join_obj)
}
