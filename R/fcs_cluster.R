#' @title  
#'   Cluster Flow Cytometry Data with Multiple Algorithms
#'
#' @description  
#'   Applies clustering to a joined flow cytometry object using one of several  
#'   algorithms: Leiden, Louvain, FlowSOM, or PhenoGraph. Can operate on raw expression data or PCA coordinates. Automatically  
#'   uses batch‐corrected data if available.
#'
#' @param fcs_join_obj  
#'   A list returned by FCSimple::fcs_join(), containing at least elements  
#'   `data` (numeric matrix) and, if `use_rep = "pca"`, `pca$pca_data`. If you’ve  
#'   already run FCSimple::fcs_batch_correction(), this function will detect  
#'   and use `fcs_join_obj$batch_correction$data`.
#'
#' @param use_rep  
#'   Character; representation to cluster. `"data"` (default) uses raw expression  
#'   values, `"pca"` uses principal‐component coordinates.
#'
#' @param language  
#'   Character; runtime environment for Python-backed clustering when applicable. `"R"` (default) runs  
#'   algorithms via R packages, `"Python"` calls out to bundled Python scripts for Leiden/Louvain.
#'
#' @param algorithm  
#'   Character; clustering algorithm to apply. One of `"leiden"`, `"louvain"`,  
#'   `"flowsom"`, or `"phenograph"`.
#'
#' @param leiden_louvain_resolution  
#'   Numeric; resolution parameter for Leiden and Louvain clustering (default 1).
#'
#' @param flowsom_nClus  
#'   Integer; number of metaclusters for FlowSOM (default 15).
#'
#' @param phenograph_k  
#'   Integer; number of nearest neighbors for PhenoGraph (default 30).
#'
#' @param adjacency_knn  
#'   Integer; k for nearest‐neighbor graph in Leiden/Louvain (default 30).
#'
#'

#' @param search_only  
#'   Logical; if `TRUE`, stops after neighbor search and returns updated  
#'   `fcs_join_obj` with either `adjacency_matrix` or `search` element  
#'   (default `FALSE`).
#'
#' @param num_cores  
#'   Integer; number of CPU cores for parallel neighbor‐search (default  
#'   `ceiling(parallel::detectCores()/2)`).
#'
#' @param output_as  
#'   Character; format of neighbor‐search output: `"adjacency"` (default) or  
#'   `"search"`.
#'
#' @details  
#'   * If batch correction is found, clustering runs on `fcs_join_obj$batch_correction$data`.  
#'   * Leiden/Louvain builds a sparse k‐NN graph and calls igraph::cluster_*().  
#'   * FlowSOM transforms data into a flowFrame and calls FlowSOM::FlowSOM().  
#'   * PhenoGraph invokes Rphenograph::Rphenograph() and extracts membership.  
#'   * After clustering, labels and algorithm settings are stored under  
#'     `fcs_join_obj[[algorithm]]`, and `object_history` is appended.
#'
#' @return  
#'   The original `fcs_join_obj` augmented with:  
#'   - For `"search_only"`: an `adjacency_matrix` or `search` element.  
#'   - For clustering algorithms: an element named by `algorithm` containing  
#'     `clusters` (factor/integer vector) and `settings` (list).  
#'   - Updated `object_history` with timestamped entry.
#'
#' @examples
#' \dontrun{
#'   ff_list <- list(flowFrame1, flowFrame2)
#'   joined <- FCSimple::fcs_join(ff_list)
#'
#'   # Leiden clustering on raw data
#'   out1 <- FCSimple::fcs_cluster(
#'     joined,
#'     use_rep = "data",
#'     algorithm = "leiden",
#'     adjacency_knn = 50,
#'     leiden_louvain_resolution = 0.5
#'   )
#'
#'   # FlowSOM clustering on PCA coordinates
#'   pca_obj <- FCSimple::fcs_pca(joined)
#'   out2 <- FCSimple::fcs_cluster(
#'     pca_obj,
#'     use_rep = "pca",
#'     algorithm = "flowsom",
#'     flowsom_nClus = 20
#'   )
#' }
#'
#' @seealso  
#'   FCSimple::fcs_join, FCSimple::fcs_batch_correction, FCSimple::fcs_pca  
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom RANN nn2
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom igraph graph.adjacency cluster_leiden cluster_louvain
#' @importFrom FlowSOM FlowSOM GetMetaclusters
#' @importFrom flowCore flowFrame
#' @importFrom Rphenograph Rphenograph
#' @export
fcs_cluster <- function(fcs_join_obj,
                        use_rep = "data", # 'data' or 'pca'
                        language = "R",
                        algorithm = c("leiden","flowsom","louvain","phenograph"),
                        leiden_louvain_resolution = 1,
                        flowsom_nClus = 15,
                        phenograph_k = 30,
                        adjacency_knn = 30,
                        search_only = FALSE,
                        num_cores = ceiling(parallel::detectCores()/2),
                        output_as = "adjacency")
{
  use_rep <- tolower(use_rep)
  if('batch_correction' %in% names(fcs_join_obj)) {
    cl_data <- fcs_join_obj[['batch_correction']][['data']]
    print("batch_correction found in fcs_join_obj list. Using fcs_join_obj[['batch_correction']][['data']] for clustering.")
  } else {
    print("batch_correction not found in fcs_join_obj list. Proceeding according to 'use_rep'.")
    if(!use_rep %in% c("data","pca")) {
      stop("'use_rep' indicates what representation of the data will be passed to the algorithm. Use 'data' or 'pca'. 'pca' requires that FCSimple::fcs_pca was already run on the data.")
    } else {
      if(use_rep=="data") {
        cl_data <- fcs_join_obj[["data"]]
        print("Using fcs_join_obj[['data']] for clustering.")
      } else if(use_rep=="pca") {
        cl_data <- fcs_join_obj[["pca"]][['pca_data']]
        print("Using fcs_join_obj[['pca']][['pca_data']] for clustering.")
      }
    }
  }
  capture_dir <- system.file(package = "FCSimple")
  if(any(length(language)!=1, !tolower(language) %in% c("r","python"))) {
    stop("error in argument 'language': use 'R' or 'Python'")
  }
  
  if(tolower(algorithm) %in% c("leiden","louvain")) {
    require(Matrix)
    if("adjacency_matrix" %in% names(fcs_join_obj)) {
      print("Adjacency matrix found, skipping nearest neighbor step.")
      sm <- fcs_join_obj[["adjacency_matrix"]]
    } else {
      require(RANN)
      require(parallel)
      require(future)
      require(future.apply)

        if(num_cores > parallel::detectCores()) {
          warning(paste0(num_cores," specified but only ",parallel::detectCores()," available. Proceeding with max available cores."))
          num_core <- parallel::detectCores()
        } else if(num_cores == 0) {
          num_core <- parallel::detectCores()
        } else {
          num_core <- num_cores
        }
        
        # Warn if too many cores specified - diminishing returns and overhead beyond ~16
        if(num_core > 16) {
          warning(paste0("Using ", num_core, " cores may cause excessive overhead. Consider using 8-16 cores for optimal performance."))
        }
        
        if(num_core > 1) {
          print(paste0("searching with ", num_core, " cores"))
        } else {
          print(paste0("searching with ", num_core, " core"))
        }
        options(future.globals.maxSize= Inf)
        future::plan("multisession", workers = num_core)

        num_neighbors <- adjacency_knn + 1
        sub_data <- vector(mode = "list", length = num_core)
        split_sums <- round(seq(from = 1, to = nrow(cl_data), length.out = num_core + 1), 0)
        
        for(i in seq_along(sub_data)) {
          if(i == 1) {
            sub_data[[i]] <- cl_data[1:(split_sums[i + 1]), , drop = FALSE]
          } else {
            sub_data[[i]] <- cl_data[(split_sums[i] + 1):(split_sums[i + 1]), , drop = FALSE]
          }
        }
        
        search_out <- NULL
        parallel_error <- NULL
        
        search_out <- tryCatch({
          future.apply::future_lapply(sub_data, FUN = function(x) {
            RANN::nn2(data = cl_data, query = x, k = num_neighbors,
                      treetype = "kd", searchtype = "standard")
          })
        }, error = function(e) {
          parallel_error <<- paste("multisession future_lapply failed:", conditionMessage(e))
          warning(parallel_error)
          NULL
        })
        
        if (is.null(search_out)) {
          # fallback to sequential with captured message in parallel_error
          if (is.null(parallel_error)) {
            parallel_error <- "multisession future_lapply returned NULL; retrying sequentially."
            warning(parallel_error)
          }
          future::plan(future::sequential)
          search_out <- lapply(sub_data, FUN = function(x) {
            RANN::nn2(data = cl_data, query = x, k = num_neighbors,
                      treetype = "kd", searchtype = "standard")
          })
        }
        
        # Reset to sequential only after all work is complete
        future::plan(future::sequential)
        gc()

        for(i in 1:length(search_out)) {
          if(i==1) {
            search_id <- search_out[[i]][[1]]
          } else {
            search_id <- rbind(search_id,search_out[[i]][[1]])
          }
        }

      nn_idx <- search_id
      if(output_as=="adjacency") {
        i_input <- rep(1:nrow(nn_idx),times=num_neighbors-1)
        j_input <- as.vector(nn_idx[,2:num_neighbors])
        sm <- Matrix::sparseMatrix(i=i_input,j=j_input,dims=c(nrow(nn_idx),nrow(nn_idx)))
        fcs_join_obj[["adjacency_matrix"]] <- sm
      } else if(output_as=="search"){
        fcs_join_obj[["search"]] <- nn_idx
      }
      if(search_only) {
        if(!'object_history' %in% names(fcs_join_obj)) {
          print("Consider running FCSimple::fcs_update() on the object.")
        }
        try(expr = fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0("nn search on ",use_rep,": ",Sys.time())), silent = TRUE)
        return(fcs_join_obj)
      }
    }
    if(tolower(language)=="python") {
      capture_dir <- system.file(package = "FCSimple") # points to package location
      Matrix::writeMM(obj = sm, file = paste0(capture_dir,"/temp_files/__python_cl_input__.mtx"))

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
      if(algorithm=="leiden") {
        fcs_join_obj[["leiden"]] <- list(clusters = cluster_numbers,
                                         settings = list(method = 'la.RBConfigurationVertexPartition',
                                                         resolution_parameter = leiden_louvain_resolution,
                                                         seed = 123, language = language))
      } else if(tolower(algorithm)=="louvain") {
        fcs_join_obj[["louvain"]] <- list(clusters = cluster_numbers,
                                          settings = list(function_call = "graph_obj.community_multilevel()",
                                                          language = language)) # left off here
      }
      if(!'object_history' %in% names(fcs_join_obj)) {
        print("Consider running FCSimple::fcs_update() on the object.")
      }
      try(expr = fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(tolower(algorithm)," on ",use_rep,": ",Sys.time())), silent = TRUE)
      return(fcs_join_obj)
    } else if(tolower(language)=="r") {
      require(igraph)
      # Symmetrize the sparse matrix for undirected graph
      sm <- sm + Matrix::t(sm)
      G <- igraph::graph.adjacency(adjmatrix = sm, mode = "undirected")
      if(tolower(algorithm)=="leiden") {
        set.seed(123)
        leid <- igraph::cluster_leiden(graph = G, objective_function = "modularity", weights = NA, resolution_parameter = leiden_louvain_resolution)
        fcs_join_obj[["leiden"]] <- list(clusters = factor(leid$membership),
                                         settings = list(resolution_parameter = leiden_louvain_resolution,
                                                         weights = NA, seed = 123, language = language))
      } else if(tolower(algorithm)=="louvain") {
        set.seed(123)
        louv <- igraph::cluster_louvain(graph = G, weights = NA, resolution = leiden_louvain_resolution)
        fcs_join_obj[["louvain"]] <- list(clusters = factor(louv$membership),
                                          settings = list(resolution = leiden_louvain_resolution,
                                                          weights = NA, seed = 123, language = language))
      }
    }
  } else {
    if(tolower(algorithm)=="flowsom") {
      require(FlowSOM)
      require(flowCore)
      som_fcs <- new(Class = "flowFrame", exprs = cl_data)
      som <- FlowSOM::FlowSOM(input = som_fcs, compensate = FALSE, transform = FALSE, silent = TRUE, nClus = flowsom_nClus)
      som_meta <- FlowSOM::GetMetaclusters(fsom = som)
      fcs_join_obj[["flowsom"]] <- list(clusters = som_meta,
                                        settings = list(compensate = FALSE, transform = FALSE,
                                                        silent = TRUE, nClus = flowsom_nClus))
      } else if(tolower(algorithm)=="phenograph") {
        require(Rphenograph)
        phenog <- Rphenograph::Rphenograph(data = cl_data, k = phenograph_k)
        phcl <- membership(phenog[[2]])
        fcs_join_obj[["phenograph"]] <- list(clusters = phcl, settings = list(k = phenograph_k))
      }
  }
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  }
  try(expr = fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(tolower(algorithm)," on ",use_rep,": ",Sys.time())), silent = TRUE)
  return(fcs_join_obj)
}

