fcs_reduce_dimensions <- function(fcs_join_obj,
                                  use_rep = "data",
                                  algorithm = c("tsne","umap"),
                                  language = c("R","Python"),
                                  umap_nn = 30,
                                  umap_min_dist = 0.1,
                                  tsne_perplexity = 30,
                                  nthread = ceiling(parallel::detectCores()/2))
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
  if(length(algorithm)!=1) {
    stop("error in argument 'algorithm': use either 'tsne' or 'umap'")
  }
  if(length(language)!=1) {
    stop("error in argument 'language': use 'R' or 'Python'")
  }
  if(tolower(algorithm)=="umap") {
    if(tolower(language)=="r") {
      require(uwot)
      require(parallel)
      set.seed(123)
      map <- uwot::umap(X = red_data, n_neighbors = round(umap_nn,0),
                        init = "spca", min_dist = umap_min_dist,
                        n_threads = ceiling(detectCores()/2), verbose = TRUE)
      colnames(map) <- c("UMAP1","UMAP2")
    } else if(tolower(language)=="python") {
      capture_dir <- system.file(package = "FCSimple")
      write.csv(red_data, file = paste0(capture_dir,"/temp_files/__python_umap_input__.csv"), row.names = FALSE)
      system(command = paste0("python ",paste0(capture_dir,"/python/run_umap.py")," ",
                              paste0(capture_dir,"/temp_files/__python_umap_input__.csv")," ",
                              capture_dir,"/temp_files ",round(umap_nn,0)," ",umap_min_dist))
      map <- read.csv(paste0(capture_dir,"/temp_files/__tmp_umap__.csv"), check.names = FALSE)
      temp_files <- list.files(path = paste0(system.file(package = "FCSimple"),"/temp_files/"), full.names = TRUE, recursive = TRUE)
      if(length(temp_files)!=0) {
        file.remove(temp_files)
      }
    } else {
      stop("error in argument 'language': use 'R' or 'Python'")
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      require(Rtsne)
      require(parallel)
      map_calculate <- Rtsne::Rtsne(X = red_data, check_duplicates = FALSE,
                                    max_iter = 2000, normalize = FALSE, perplexity = round(tsne_perplexity,0),
                                    stop_lying_iter = 700, mom_switch_iter = 700,
                                    eta = round(nrow(red_data)/12),
                                    num_threads = nthread)
      map <- map_calculate[["Y"]]
      colnames(map) <- c("tSNE1","tSNE2")
    } else if(tolower(language)=="python") {
      require(parallel)
      capture_dir <- system.file(package = "FCSimple")
      write.csv(red_data, file = paste0(capture_dir,"/temp_files/__python_tsne_input__.csv"), row.names = FALSE)
      system(command = paste0("python ",paste0(capture_dir,"/python/run_tsne.py")," ",
                              paste0(capture_dir,"/temp_files/__python_tsne_input__.csv")," ",
                              capture_dir,"/temp_files"," ",floor(parallel::detectCores()/2)," ",round(tsne_perplexity,0)))
      map <- read.csv(paste0(capture_dir,"/temp_files/__tmp_tsne__.csv"), check.names = FALSE)
      temp_files <- list.files(path = paste0(system.file(package = "FCSimple"),"/temp_files/"), full.names = TRUE, recursive = TRUE)
      if(length(temp_files)!=0) {
        file.remove(temp_files)
      }
    } else {
      stop("error in argument 'language': use 'R' or 'Python'")
    }
  }
  coordinates_list <- map
  if(tolower(algorithm)=="umap") {
    if(tolower(language)=="r") {
      settings_list <- list(use_rep = use_rep, language = "R", init = "spca", min_dist = 0.1,
                            n_threads = ceiling(detectCores()/2), num_neighbors = round(umap_nn,0),
                            min_dist = umap_min_dist, verbose = TRUE)
    } else if(tolower(language)=="python") {
      settings_list <- list(use_rep = use_rep, language = "Python", init = 'spectral', low_memory = 'True',
                            random_state = 123, num_neighbors = round(umap_nn,0),
                            min_dist = umap_min_dist, transform_seed = 123, verbose = 'True')
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      settings_list <- list(use_rep = use_rep, language = "R", check_duplicates = FALSE, max_iter = 2000,
                            normalize = FALSE, stop_lying_iter = 700, mom_switch_iter = 700,
                            eta = round(nrow(map_input)/12), perplexity = round(tsne_perplexity,0),
                            num_threads = ceiling(detectCores()/2))
    }
    if(tolower(language)=="python") {
      settings_list <- list(use_rep = use_rep, language = "Python", perplexity = 30,
                            metric = "euclidean", random_state = 123, verbose = "True",
                            perplexity = round(tsne_perplexity,0),
                            num_threads = ceiling(detectCores()/2))
    }
  }
  fcs_join_obj[[length(fcs_join_obj)+1]] <- list(coordinates = coordinates_list,
                                                 settings = settings_list)
  names(fcs_join_obj)[length(fcs_join_obj)] <- ifelse(tolower(algorithm)=="umap","umap","tsne")
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  }
  try(expr = fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']], paste0(tolower(algorithm)," on ",use_rep,": ",Sys.time())), silent = TRUE)
  return(fcs_join_obj)
}
