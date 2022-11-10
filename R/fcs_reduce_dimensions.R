fcs_reduce_dimensions <- function(fcs_join_obj,
                                  algorithm = c("tsne","umap"),
                                  language = c("R","Python"))
{
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
      map <- uwot::umap(X = fcs_join_obj[["data"]], n_neighbors = 30, init = "spca", min_dist = 0.1,
                        n_threads = ceiling(detectCores()/2), verbose = TRUE)
      colnames(map) <- c("UMAP1","UMAP2")
    } else if(tolower(language)=="python") {
      capture_dir <- system.file(package = "FCSimple")
      write.csv(fcs_join_obj[["data"]], file = paste0(capture_dir,"/temp_files/__python_umap_input__.csv"), row.names = FALSE)
      system(command = paste0("python ",paste0(capture_dir,"/python/run_umap.py")," ",paste0(capture_dir,"/temp_files/__python_umap_input__.csv")," ",capture_dir,"/temp_files"))
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
      map_calculate <- Rtsne::Rtsne(X = fcs_join_obj[["data"]], check_duplicates = FALSE, max_iter = 2000, normalize = FALSE,
                                    stop_lying_iter = 700, mom_switch_iter = 700,
                                    eta = round(nrow(map_input)/12),
                                    num_threads = ceiling(detectCores()/2))
      map <- map_calculate[["Y"]]
      colnames(map) <- c("tSNE1","tSNE2")
    } else if(tolower(language)=="python") {
      require(parallel)
      capture_dir <- system.file(package = "FCSimple")
      write.csv(fcs_join_obj[["data"]], file = paste0(capture_dir,"/temp_files/__python_tsne_input__.csv"), row.names = FALSE)
      system(command = paste0("python ",paste0(capture_dir,"/python/run_tsne.py")," ",paste0(capture_dir,"/temp_files/__python_tsne_input__.csv")," ",capture_dir,"/temp_files"," ",floor(parallel::detectCores()/2)))
      map <- read.csv(paste0(capture_dir,"/temp_files/__tmp_tsne__.csv"), check.names = FALSE)
      temp_files <- list.files(path = paste0(system.file(package = "FCSimple"),"/temp_files/"), full.names = TRUE, recursive = TRUE)
      if(length(temp_files)!=0) {
        file.remove(temp_files)
      }
    }
  }
  coordinates_list <- map
  if(tolower(algorithm)=="umap") {
    if(tolower(language)=="r") {
      settings_list <- list(language = "R", n_neighbors = 30, init = "spca", min_dist = 0.1,
                            n_threads = ceiling(detectCores()/2), verbose = TRUE)
    } else if(tolower(language)=="python") {
      settings_list <- list(language = "Python", n_neighbors = 30, init = 'spectral',
                            min_dist = 0.1, low_memory = 'True', random_state = 123,
                            transform_seed = 123, verbose = 'True')
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      settings_list <- list(language = "R", check_duplicates = FALSE, max_iter = 2000,
                            normalize = FALSE, stop_lying_iter = 700, mom_switch_iter = 700,
                            eta = round(nrow(map_input)/12), num_threads = ceiling(detectCores()/2))
    }
    if(tolower(language)=="python") {
      settings_list <- list(language = "Python", perplexity = 30, metric = "euclidean",
                            random_state = 123, verbose = "True", num_threads = ceiling(detectCores()/2))
    }
  }
  fcs_join_obj[[length(fcs_join_obj)+1]] <- list(coordinates = coordinates_list,
                                                 settings = settings_list)
  names(fcs_join_obj)[length(fcs_join_obj)] <- ifelse(tolower(algorithm)=="umap","umap","tsne")
  return(fcs_join_obj)
}
