#' @title Batch Correction of Joined Flow Cytometry Data
#'
#' @description
#' Applies batch‐effect correction to a joined flow cytometry object using
#' either cyCombine or Harmony. The function can operate on raw data or PCA
#' coordinates and stores corrected values along with metadata and history.
#'
#' @param fcs_join_obj A list returned by FCSimple::fcs_join(), containing at
#'   minimum elements "data", "run_date", and "source". If `use_rep = "pca"`,
#'   it must also contain `fcs_join_obj$pca$pca_data`.
#'
#' @param use_rep Character; either "data" (default) to correct raw expression
#'   values, or "pca" to correct principal‐component coordinates.
#'
#' @param correction_method Character vector; one of "cyCombine" (default) or
#'   "harmony", selecting the algorithm for batch correction.
#'
#' @param correction_markers Character vector of channel names to correct
#'   (default "all"). Only used with `correction_method = "cyCombine"`.
#'
#' @param batch_source_regex Regular expression (string) to extract batch IDs
#'   from the "source" element. Default `"[0-9]+\\-[A-Za-z]+\\-[0-9]+"`.
#'
#' @param cyCombine_SOMx Integer; number of grid columns for the SOM in
#'   cyCombine (default 8).
#'
#' @param cyCombine_SOMy Integer; number of grid rows for the SOM in
#'   cyCombine (default 8).
#'
#' @param cyCombine_detect_effects Logical; if `TRUE`, runs cyCombine’s
#'   detect_batch_effects before and after normalization (default FALSE).
#'
#' @param harmony_cores Integer; number of CPU cores for Harmony (default 1).
#'
#' @param harmony_iterations Integer; maximum number of Harmony iterations
#'   (default 10).
#'
#' @param harmony_covars Character vector of metadata columns in
#'   `fcs_join_obj` to use as covariates in Harmony (default `c("batch","sample")`).
#'
#' @param harmony_lambda Numeric; ridge‐penalty parameter for Harmony to
#'   prevent overcorrection (default 1).
#'
#' @param harmony_sample_element Name of the element in `fcs_join_obj` to treat
#'   as the sample identifier for Harmony (default "source").
#'
#' @details
#' - For `cyCombine`, raw data are scaled, SOM clusters created, and per‐batch
#'   normalization performed. Corrected expression values are stored under
#'   `fcs_join_obj$batch_correction$data` along with a history list.
#' - For `harmony`, PCA or raw data matrix is passed to `RunHarmony()`, and
#'   corrected embeddings plus metadata are returned in
#'   `fcs_join_obj$batch_correction`.
#'
#' @return
#' The input `fcs_join_obj` augmented with:
#' - `batch_correction`: a list containing corrected data, method, and metadata.
#' - `object_history`: appended entry recording the batch correction event.
#'
#' @examples
#' \dontrun{
#' ff_list <- list(flowFrame1, flowFrame2)
#' joined <- FCSimple::fcs_join(ff_list, run_date = c("2023-01-01","2023-01-02"))
#'
#' corrected1 <- FCSimple::fcs_batch_correction(
#'   joined,
#'   use_rep = "data",
#'   correction_method = "cyCombine",
#'   cyCombine_detect_effects = TRUE
#' )
#'
#' pca_obj <- FCSimple::fcs_pca(joined)
#' corrected2 <- FCSimple::fcs_batch_correction(
#'   pca_obj,
#'   use_rep = "pca",
#'   correction_method = "harmony",
#'   harmony_cores = 4,
#'   harmony_covars = c("batch", "sample")
#' )
#' }
#'
#' @seealso
#' FCSimple::fcs_join, FCSimple::fcs_pca, cyCombine::normalize, harmony::RunHarmony
#'
#' @importFrom cyCombine normalize get_markers detect_batch_effects
#' @importFrom harmony RunHarmony
#' @importFrom stringr str_extract
#' @export
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
