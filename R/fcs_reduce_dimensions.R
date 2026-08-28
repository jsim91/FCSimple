#' @title 2D Dimensionality Reduction for Flow Cytometry Data
#'
#' @description
#'   Performs UMAP or t‐SNE on a flow cytometry analysis object. By default,
#'   reduces the (batch‐corrected) expression matrix to two dimensions
#'   using either the R implementation (uwot or Rtsne) or an external
#'   Python script. The result is stored in your object under “umap” or “tsne”.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join() (and optionally
#'   augmented by FCSimple::fcs_batch_correction or FCSimple::fcs_pca). Must
#'   contain at least one of:
#'   - `fcs_join_obj$data` (raw or transformed expression matrix), or
#'   - `fcs_join_obj$pca$pca_data` (when `use_rep = "pca"`), or
#'   - `fcs_join_obj$batch_correction$data` (automatically used if present).
#'
#' @param use_rep
#'   Character; which representation to reduce.
#'   - `"data"` (default): use `fcs_join_obj$data` or batch‐corrected data if present.
#'   - `"pca"`: use `fcs_join_obj$pca$pca_data` (requires prior call to FCSimple::fcs_pca).
#'
#' @param algorithm
#'   Character; which algorithm to run.
#'   - `"tsne"`: t‐distributed stochastic neighbor embedding.
#'   - `"umap"`: uniform manifold approximation and projection.
#'   Default is `c("tsne","umap")` (selects first).
#'
#' @param language
#'   Character; runtime environment for the chosen algorithm.
#'   - `"R"` (default): calls uwot::umap or Rtsne::Rtsne.
#'   - `"Python"`: writes a CSV, invokes the package’s Python script, and reads
#'     back results. For t‐SNE, the Python backend runs the opt‐SNE
#'     implementation (omiq‐ai Multicore‐opt‐SNE) in automated‐parameter mode,
#'     or falls back to openTSNE automatically when opt‐SNE is not importable.
#'
#' @param umap_nn
#'   Numeric; number of neighbors for UMAP (default 30).
#'
#' @param umap_min_dist
#'   Numeric; minimum distance parameter for UMAP (default 0.1).
#'
#' @param tsne_perplexity
#'   Numeric; perplexity parameter for t‐SNE (default 30).
#'
#' @param n_components
#'   Integer; number of output dimensions for the reduction. Default `2`
#'   (standard 2D embedding). Set `3` for a 3D embedding suitable for
#'   interactive visualization with FCSimple::fcs_plot_reduction_3d().
#'
#' @param num_cores
#'   Integer; number of CPU threads for parallel computation (default `ceiling(parallel::detectCores()/2)`).
#'
#' @param seed
#'   Integer or NA; random seed for reproducibility (default NA).
#'   When NA, Python implementations will not set a seed, allowing for faster
#'   multi-threaded computation. When set to an integer, both R and Python
#'   implementations will use this seed for reproducible results.
#'
#' @details
#'   1. If `fcs_join_obj$batch_correction$data` exists, that matrix is used
#'      regardless of `use_rep`. Otherwise, `use_rep` selects raw data or PCA.
#'   2. For UMAP:
#'      - R: calls `uwot::umap()` with `umap_nn` and `umap_min_dist`; uses `seed` if provided.
#'      - Python: writes data to `inst/python`, runs `run_umap.py`, cleans temp files.
#'   3. For t‐SNE:
#'      - R: calls `Rtsne::Rtsne()` with `tsne_perplexity`, `num_cores`, and
#'        fixed settings; an advisory message reports opt‐SNE readiness.
#'      - Python: runs opt‐SNE (auto-parameterized t‐SNE) via `run_optsne.py`
#'        when available, otherwise falls back to openTSNE via `run_tsne.py`.
#'   4. The resulting 2‐column matrix is stored as
#'      `fcs_join_obj$umap$coordinates` or `$tsne$coordinates`, and the
#'      parameters used are recorded under
#'      `fcs_join_obj$<algorithm>$settings`.
#'   5. An entry is appended to `fcs_join_obj$object_history`:
#'      `<algorithm> on <use_rep>: <timestamp>`.
#'
#' @return
#'   The input `fcs_join_obj`, with a new element named by the
#'   lower‐case `algorithm`:
#'   - `$<algorithm>$coordinates`: numeric matrix (cells × 2).
#'   - `$<algorithm>$settings`: list of parameters passed, including a
#'     `features` element naming the features used for the reduction.
#'   - `object_history` updated with the reduction event.
#'
#' @examples
#' \dontrun{
#'   # Basic UMAP on raw data
#'   joined <- FCSimple::fcs_join(files)
#'   out_umap <- FCSimple::fcs_reduce_dimensions(
#'     joined,
#'     algorithm = "umap",
#'     language  = "R"
#'   )
#'
#'   # t-SNE (opt-SNE) using PCA coordinates and Python backend
#'   pca_obj <- FCSimple::fcs_pca(joined)
#'   out_tsne <- FCSimple::fcs_reduce_dimensions(
#'     pca_obj,
#'     use_rep    = "pca",
#'     algorithm  = "tsne",
#'     language   = "Python",
#'     tsne_perplexity = 50
#'   )
#' }
#'
#' @seealso
#'   uwot::umap, Rtsne::Rtsne, FCSimple::fcs_pca, FCSimple::fcs_batch_correction.
#'   opt‐SNE: https://github.com/omiq-ai/Multicore-opt-SNE
#'
#' @importFrom uwot umap
#' @importFrom Rtsne Rtsne
#' @importFrom parallel detectCores
#' @export
fcs_reduce_dimensions <- function(fcs_join_obj,
                                  use_rep = "data",
                                  algorithm = c("tsne","umap"),
                                  language = "R",
                                  umap_nn = 30,
                                  umap_min_dist = 0.1,
                                  tsne_perplexity = 30,
                                  n_components = 2,
                                  num_cores = ceiling(parallel::detectCores()/2),
                                  seed = NA)
{
  use_rep <- tolower(use_rep)
  if('batch_correction' %in% names(fcs_join_obj)) {
    red_data <- fcs_join_obj[['batch_correction']][['data']]
    print("batch_correction found in fcs_join_obj list. Using fcs_join_obj[['batch_correction']][['data']] for dimension reduction.")
  } else {
    print("batch_correction not found in fcs_join_obj list. Proceeding according to 'use_rep'.")
    if(!use_rep %in% c("data","pca")) {
      stop("'use_rep' indicates what representation of the data will be passed to the algorithm. Use 'data' or 'pca'. 'pca' requires that FCSimple::fcs_pca was already run on the data.")
    } else {
      if(use_rep=="data") {
        red_data <- fcs_join_obj[["data"]]
        print("Using fcs_join_obj[['data']] for dimension reduction.")
      } else if(use_rep=="pca") {
        red_data <- fcs_join_obj[["pca"]][['pca_data']]
        print("Using fcs_join_obj[['pca']][['pca_data']] for dimension reduction.")
      }
    }
  }
  red_features <- .fcs_features(red_data)
  if(length(algorithm)!=1) {
    stop("error in argument 'algorithm': use either 'tsne' or 'umap'")
  }
  if(length(language)!=1) {
    stop("error in argument 'language': use 'R' or 'Python'")
  }
  if(tolower(algorithm)=="umap") {
    if(tolower(language)=="r") {
      if (!require(uwot, quietly = TRUE)) stop("Package 'uwot' is required but could not be loaded.")
      if (!require(parallel, quietly = TRUE)) stop("Package 'parallel' is required but could not be loaded.")
      if(!is.na(seed)) {
        set.seed(seed)
      }
      map <- uwot::umap(X = red_data, n_neighbors = round(umap_nn,0),
                        n_components = n_components,
                        init = "spca", min_dist = umap_min_dist,
                        n_threads = num_cores, verbose = TRUE)
      colnames(map) <- paste0("UMAP", seq_len(n_components))
    } else if(tolower(language)=="python") {
      temp_dir <- .fcs_temp_dir()
      on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
      umap_script <- system.file("python", "run_umap.py", package = "FCSimple")
      umap_in <- file.path(temp_dir, "__python_umap_input__.csv")
      umap_out <- file.path(temp_dir, "__tmp_umap__.csv")
      data.table::fwrite(data.table::as.data.table(red_data), file = umap_in,
                         nThread = parallel::detectCores(), row.names = FALSE)
      umap_exit <- system(command = paste("python", shQuote(umap_script), shQuote(umap_in), shQuote(temp_dir),
                                          round(umap_nn,0), umap_min_dist, num_cores, n_components))
      if(!identical(umap_exit, 0L)) stop("Python UMAP failed with exit code ", umap_exit)
      if(!file.exists(umap_out)) stop("Python UMAP did not produce an output file.")
      map <- read.csv(umap_out, check.names = FALSE)
    } else {
      stop("error in argument 'language': use 'R' or 'Python'")
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      if (!require(Rtsne, quietly = TRUE)) stop("Package 'Rtsne' is required but could not be loaded.")
      if (!require(parallel, quietly = TRUE)) stop("Package 'parallel' is required but could not be loaded.")
      .emit_optsne_advisory()
      if(!is.na(seed)) {
        set.seed(seed)
      }
      map_calculate <- Rtsne::Rtsne(X = red_data, check_duplicates = FALSE,
                                    dims = n_components,
                                    max_iter = 2000, normalize = FALSE, perplexity = round(tsne_perplexity,0),
                                    stop_lying_iter = 700, mom_switch_iter = 700,
                                    eta = round(nrow(red_data)/12),
                                    num_threads = num_cores)
      map <- map_calculate[["Y"]]
      colnames(map) <- paste0("tSNE", seq_len(n_components))
    } else if(tolower(language)=="python") {
      temp_dir <- .fcs_temp_dir()
      on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
      tsne_in <- file.path(temp_dir, "__python_tsne_input__.csv")
      tsne_out <- file.path(temp_dir, "__tmp_tsne__.csv")
      write.csv(red_data, file = tsne_in, row.names = FALSE)
      seed_arg <- ifelse(is.na(seed), "NA", as.character(seed))

      # Prefer opt-SNE; automatically fall back to openTSNE when the bundled
      # MulticoreTSNE extension is not importable.
      tsne_python_backend <- if (isTRUE(.check_optsne_prereqs()$optsne_importable)) {
        "opt-SNE"
      } else {
        message("opt-SNE prerequisites not met; falling back to openTSNE for t-SNE.")
        "openTSNE"
      }
      tsne_script <- if (tsne_python_backend == "opt-SNE") "run_optsne.py" else "run_tsne.py"
      tsne_script_path <- system.file("python", tsne_script, package = "FCSimple")

      tsne_exit <- system(command = paste("python", shQuote(tsne_script_path), shQuote(tsne_in), shQuote(temp_dir),
                                          floor(parallel::detectCores()/2), round(tsne_perplexity,0), seed_arg, n_components))
      if(!identical(tsne_exit, 0L)) stop("Python t-SNE failed with exit code ", tsne_exit)
      if(!file.exists(tsne_out)) stop("Python t-SNE did not produce an output file.")
      map <- read.csv(tsne_out, check.names = FALSE)
    } else {
      stop("error in argument 'language': use 'R' or 'Python'")
    }
  }
  coordinates_list <- map
  if(tolower(algorithm)=="umap") {
    if(tolower(language)=="r") {
      settings_list <- list(use_rep = use_rep, language = "R", init = "spca",
                            n_components = n_components,
                            n_threads = ceiling(detectCores()/2), num_neighbors = round(umap_nn,0),
                            min_dist = umap_min_dist, verbose = TRUE, seed = seed)
    } else if(tolower(language)=="python") {
      if(is.na(seed)) {
        settings_list <- list(use_rep = use_rep, language = "Python", init = 'spectral', low_memory = 'True',
                              n_components = n_components,
                              num_neighbors = round(umap_nn,0),
                              min_dist = umap_min_dist, n_jobs = num_cores, verbose = 'True', seed = NA)
      } else {
        settings_list <- list(use_rep = use_rep, language = "Python", init = 'spectral', low_memory = 'True',
                              n_components = n_components,
                              random_state = seed, num_neighbors = round(umap_nn,0),
                              min_dist = umap_min_dist, transform_seed = seed, n_jobs = 1, verbose = 'True', seed = seed)
      }
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      settings_list <- list(use_rep = use_rep, language = "R", check_duplicates = FALSE, max_iter = 2000,
                            dims = n_components,
                            normalize = FALSE, stop_lying_iter = 700, mom_switch_iter = 700,
                            eta = round(nrow(red_data)/12), perplexity = round(tsne_perplexity,0),
                            num_threads = ceiling(detectCores()/2), seed = seed)
    }
    if(tolower(language)=="python") {
      if (tsne_python_backend == "opt-SNE") {
        settings_list <- list(use_rep = use_rep, language = "Python",
                              method = "opt-SNE", fallback = FALSE,
                              auto_iter = TRUE, auto_iter_end = 5000,
                              early_exaggeration = 12, angle = 0.5,
                              n_components = n_components,
                              perplexity = round(tsne_perplexity,0),
                              num_threads = ceiling(detectCores()/2), seed = seed)
      } else {
        settings_list <- list(use_rep = use_rep, language = "Python",
                              method = "openTSNE", fallback = TRUE,
                              metric = "euclidean",
                              n_components = n_components,
                              perplexity = round(tsne_perplexity,0),
                              num_threads = ceiling(detectCores()/2), seed = seed)
      }
    }
  }
  settings_list$features <- red_features
  fcs_join_obj[[length(fcs_join_obj)+1]] <- list(coordinates = coordinates_list,
                                                 settings = settings_list)
  slot_name <- paste0(ifelse(tolower(algorithm)=="umap","umap","tsne"), "_", n_components, "d")
  names(fcs_join_obj)[length(fcs_join_obj)] <- slot_name

  # Track reduction creation order so the first-created reduction drives
  # cluster colouring consistently across 2D and 3D plots.
  if (is.null(fcs_join_obj$reduction_order)) {
    fcs_join_obj$reduction_order <- slot_name
  } else {
    fcs_join_obj$reduction_order <- c(fcs_join_obj$reduction_order, slot_name)
  }
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_audit() on the object.")
  }
  try(expr = fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(tolower(algorithm)," on ",use_rep,": ",Sys.time())), silent = TRUE)
  return(fcs_join_obj)
}

# Check whether the Python opt-SNE backend is usable.
.check_optsne_prereqs <- function() {
  python_available <- nzchar(Sys.which("python"))
  cmake_available <- nzchar(Sys.which("cmake"))
  optsne_importable <- FALSE
  if (python_available) {
    probe <- system.file("python", "check_optsne.py", package = "FCSimple")
    if (nzchar(probe) && file.exists(probe)) {
      out <- tryCatch(
        suppressWarnings(system2(command = "python", args = shQuote(probe), stdout = TRUE, stderr = TRUE)),
        error = function(e) NULL
      )
      status <- attr(out, "status")
      # system2 attaches a non-NULL "status" attribute only on failure; a
      # successful run returns output with no status attribute (or status 0).
      optsne_importable <- !is.null(out) && (is.null(status) || isTRUE(status == 0))
    }
  }
  list(python_available = python_available,
       optsne_importable = optsne_importable,
       cmake_available = cmake_available)
}

# Emit an advisory message when the R t-SNE backend is chosen.
.emit_optsne_advisory <- function() {
  pr <- .check_optsne_prereqs()
  message(sprintf(
    paste0("Using Rtsne (R backend). A higher-quality alternative, opt-SNE, is available via language = 'Python'. ",
           "opt-SNE prerequisites: Python installed: %s; opt-SNE package importable: %s; cmake installed (needed to build opt-SNE): %s. ",
           "Install missing pieces with FCSimple::fcs_install_python_dependencies(install = TRUE, build_optsne = TRUE)."),
    ifelse(pr$python_available, "yes", "no"),
    ifelse(pr$optsne_importable, "yes", "no"),
    ifelse(pr$cmake_available, "yes", "no")
  ))
}
