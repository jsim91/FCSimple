fcs_cluster <- function(fcs_join_obj,
                        language = c("R","Python"),
                        algorithm = c("leiden","flowsom","louvain","phenograph"),
                        leiden_louvain_resolution = 1,
                        flowsom_nClus = 20,
                        phenograph_k = 30,
                        adjacency_knn = 30)
{
  if(any(length(language)!=1, !tolower(language) %in% c("r","python"))) {
    stop("error in argument 'language': use 'R' or 'Python'")
  }
  if(tolower(algorithm) %in% c("leiden","louvain")) {
    require(FNN)
    require(Matrix)
    if(all("nn_search" %in% names(fcs_join_obj), c("matrix","data.frame") %in% class(fcs_join_obj[["nn_search"]]))) {
      adj_search <- fcs_join_obj[["nn_search"]]
      adjacency_knn <- ncol(adj_search)
    } else {
      adj_search <- FNN::knn.index(data = fcs_join_obj[["data"]], k = adjacency_knn)
      fcs_join_obj[["nn_search"]] <- adj_search
    }
    i_input <- rep(1:nrow(adj_search),times=ncol(adj_search))
    j_input <- as.vector(adj_search)
    sm <- Matrix::sparseMatrix(i=i_input,j=j_input,dims=c(nrow(adj_search),nrow(adj_search)))
    if(tolower(language)=="python") {
      capture_dir <- system.file(package = "FCSimple")
      Matrix::writeMM(obj = sm, file = paste0(capture_dir,"/py/__python_cl_input__.mtx"))

      system(command = paste0("python ",paste0(capture_dir,"/py/run_cluster.py")," ",paste0(capture_dir,"/py/__python_cl_input__.mtx")," ",capture_dir," ",tolower(algorithm)," ",leiden_louvain_resolution))
      map <- read.csv(paste0(capture_dir,"/py/__tmp_cl__.csv"), check.names = FALSE)
      if (file.exists(paste0(capture_dir,"/py/__tmp_cl__.csv"))) {
        file.remove(paste0(capture_dir,"/py/__tmp_cl__.csv"))
      }
      cluster_numbers <- map[,1]
      if(algorithm=="leiden") {
        fcs_join_obj[["leiden"]] <- list(clusters = cluster_numbers,
                                               settings = list(method = 'la.RBConfigurationVertexPartition',
                                                               resolution_parameter = leiden_louvain_resolution,
                                                               seed = 123))
      } else if(tolower(algorithm)=="louvain") {
        fcs_join_obj[["louvain"]] <- list(clusters = cluster_numbers,
                                               function_call = "graph_obj.community_multilevel()") # left off here
      }
      return(fcs_join_obj)
    } else if(tolower(language)=="r") {
      if(tolower(algorithm)=="leiden") {
        require(leiden)
        leid <- leiden::leiden(object = sm, resolution_parameter = leiden_louvain_resolution, seed = 123)
        fcs_join_obj[["leiden"]] <- list(clusters = leid,
                                               settings = list(resolution_parameter = leiden_louvain_resolution,
                                                               seed = 123))
      } else if(tolower(algorithm)=="louvain") {
        require(igraph)
        G <- igraph::graph.adjacency(adjmatrix = sm)
        louv <- igraph::cluster_louvain(graph = G, weights = NULL, resolution = leiden_louvain_resolution)
        fcs_join_obj[["louvain"]] <- list(clusters = louv$membership,
                                               settings = list(resolution = leiden_louvain_resolution,
                                                               weights = NULL))
      }
    }
  }
  if(tolower(algorithm)=="flowsom") {
    require(FlowSOM)
    require(flowCore)
    som_fcs <- new(Class = "flowFrame", exprs = fcs_join_obj[["data"]])
    som <- FlowSOM::FlowSOM(input = som_fcs, compensate = FALSE, transform = FALSE, silent = TRUE, nClus = som_nClus)
    som_meta <- FlowSOM::GetMetaclusters(fsom = som)
    fcs_join_obj[["flowsom"]] <- list(clusters = som_meta,
                                           settings = list(compensate = FALSE, transform = FALSE,
                                                           silent = TRUE, nClus = som_nClus))
  }
  if(tolower(algorithm)=="phenograph") {
    require(Rphenograph)
    phenog <- Rphenograph::Rphenograph(data = fcs_join_obj[["data"]], k = phenograph_k)
    phcl <- membership(phenog[[2]])
    fcs_join_obj[["phenograph"]] <- list(clusters = phcl,
                                           settings = list(k = phenograph_k))
  }
  return(fcs_join_obj)
}
