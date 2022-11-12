fcs_cluster <- function(fcs_join_obj,
                        language = c("R","Python"),
                        algorithm = c("leiden","flowsom","louvain","phenograph"),
                        leiden_louvain_resolution = 1,
                        flowsom_nClus = 15,
                        phenograph_k = 30,
                        adjacency_knn = 30,
                        search_method = c("FNN","RANN"))
{
  capture_dir <- system.file(package = "FCSimple")
  if(any(length(language)!=1, !tolower(language) %in% c("r","python"))) {
    stop("error in argument 'language': use 'R' or 'Python'")
  }
  if(tolower(algorithm) %in% c("leiden","louvain")) {
    require(FNN)
    require(Matrix)
    if("adjacency_matrix" %in% names(fcs_join_obj)) {
      if(fcs_join_obj[["adjacency_matrix"]]=="ngCMatrix") {
        sm <- fcs_join_obj[["adjacency_matrix"]]
      }
    } else {
      if(tolower(search_method)=="fnn") {
        adj_search <- FNN::knn.index(data = fcs_join_obj[["data"]], k = adjacency_knn)
        i_input <- rep(1:nrow(adj_search),times=ncol(adj_search))
        j_input <- as.vector(adj_search)
        sm <- Matrix::sparseMatrix(i=i_input,j=j_input,dims=c(nrow(adj_search),nrow(adj_search)))
        fcs_join_obj[["adjacency_matrix"]] <- sm
      } else if(tolower(search_method)=="rann") {
        require(RANN)
        require(parallel)
        require(future)
        require(future.apply)

        num_core <- ceiling(detectCores()/2)
        options(future.globals.maxSize= Inf)
        future::plan(future::cluster, workers = num_core)

        num_neighbors <- adjacency_knn + 1
        sub_data <- vector(mode="list",length=num_core)
        split_sums <- round(seq(from=1,to=nrow(fcs_join_obj[["data"]]),length.out=num_core+1),0)
        for(i in 1:length(sub_data)) {
          if(i==1) {
            sub_data[[i]] <- fcs_join_obj[["data"]][1:(split_sums[i+1]),]
          } else {
            sub_data[[i]] <- fcs_join_obj[["data"]][(split_sums[i] + 1):(split_sums[i+1]),]
          }
        }
        search_out <- future_lapply(sub_data,FUN=function(x) {
          return(RANN::nn2(data=fcs_join_obj[["data"]],query=x,k=num_neighbors,treetype = "kd",searchtype = "standard"))
        })
        future::plan(future::sequential)

        for(i in 1:length(search_out)) {
          if(i==1) {
            search_id <- search_out[[i]][[1]]
          } else {
            search_id <- rbind(search_id,search_out[[i]][[1]])
          }
        }

        nn_idx <- search_id
        i_input <- rep(1:nrow(nn_idx),times=num_neighbors-1)
        j_input <- as.vector(nn_idx[,2:num_neighbors])
        sm <- Matrix::sparseMatrix(i=i_input,j=j_input,dims=c(nrow(nn_idx),nrow(nn_idx)))
        fcs_join_obj[["adjacency_matrix"]] <- sm
      } else {
        stop("error in argument 'search_method': use either 'RANN' or 'FNN'")
      }
    }
    if(tolower(language)=="python") {
      capture_dir <- system.file(package = "FCSimple") # points to package location
      Matrix::writeMM(obj = sm, file = paste0(capture_dir,"/python/__python_cl_input__.mtx"))

      system(command = paste0("python ",paste0(capture_dir,"/python/run_cluster.py")," ",
                              paste0(capture_dir,"/python/__python_cl_input__.mtx")," ",capture_dir,"/python ",tolower(algorithm)," ",leiden_louvain_resolution))
      map <- read.csv(paste0(capture_dir,"/python/__tmp_cl__.csv"), check.names = FALSE)
      if(file.exists(paste0(capture_dir,"/python/__tmp_cl__.csv"))) {
        file.remove(paste0(capture_dir,"/python/__tmp_cl__.csv"))
      }
      if(file.exists(paste0(capture_dir,"/python/__python_cl_input__.mtx"))) {
        file.remove(paste0(capture_dir,"/python/__python_cl_input__.mtx"))
      }
      cluster_numbers <- map[,1]
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
      return(fcs_join_obj)
    } else if(tolower(language)=="r") {
      G <- igraph::graph.adjacency(adjmatrix = sm, mode = "undirected")
      require(igraph)
      if(tolower(algorithm)=="leiden") {
        set.seed(123)
        leid <- igraph::cluster_leiden(graph = G, objective_function = "modularity", weights = NA, resolution_parameter = leiden_louvain_resolution, )
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
      som_fcs <- new(Class = "flowFrame", exprs = fcs_join_obj[["data"]])
      som <- FlowSOM::FlowSOM(input = som_fcs, compensate = FALSE, transform = FALSE, silent = TRUE, nClus = flowsom_nClus)
      som_meta <- FlowSOM::GetMetaclusters(fsom = som)
      fcs_join_obj[["flowsom"]] <- list(clusters = som_meta,
                                        settings = list(compensate = FALSE, transform = FALSE,
                                                        silent = TRUE, nClus = flowsom_nClus))
      } else if(tolower(algorithm)=="phenograph") {
        require(Rphenograph)
        phenog <- Rphenograph::Rphenograph(data = fcs_join_obj[["data"]], k = phenograph_k)
        phcl <- membership(phenog[[2]])
        fcs_join_obj[["phenograph"]] <- list(clusters = phcl, settings = list(k = phenograph_k))
      }
  }
  return(fcs_join_obj)
}
