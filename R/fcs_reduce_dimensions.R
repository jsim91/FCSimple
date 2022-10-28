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
      capture_dir <- system.file(package = "FCSimple") # points to package location
      write.csv(fcs_join_obj[["data"]], file = paste0(capture_dir,"/python/__python_umap_input__.csv"), row.names = FALSE)
      system(command = paste0("python ",paste0(capture_dir,"/python/run_umap.py")," ",paste0(capture_dir,"/python/__python_umap_input__.csv")," ",capture_dir,"/python"))
      map <- read.csv(paste0(capture_dir,"/python/__tmp_umap__.csv"), check.names = FALSE)
      if (file.exists(paste0(capture_dir,"/python/__tmp_umap__.csv"))) {
        file.remove(paste0(capture_dir,"/python/__tmp_umap__.csv"))
      }
      if (file.exists(paste0(capture_dir,"/python/__python_umap_input__.csv"))) {
        file.remove(paste0(capture_dir,"/python/__python_umap_input__.csv"))
      }
    } else {
      stop("error in argument 'language': use 'R' or 'Python'")
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      require(Rtsne)
      require(parallel)
      map_input <- Rtsne::normalize_input(fcs_join_obj[["data"]])
      map_calculate <- Rtsne::Rtsne(X = map_input, check_duplicates = FALSE, max_iter = 2000, normalize = FALSE,
                                    stop_lying_iter = 700, mom_switch_iter = 700,
                                    eta = round(nrow(map_input)/12),
                                    num_threads = ceiling(detectCores()/2))
      map <- map_calculate[["Y"]]
      colnames(map) <- c("tSNE1","tSNE2")
    } else if(tolower(language)=="python") {
      warning("error in 'language' + 'algorithm': tsne in python not supported. Using R for now.")
    }
    require(Rtsne)
    require(parallel)
    map_input <- Rtsne::normalize_input(fcs_join_obj[["data"]])
    map_calculate <- Rtsne::Rtsne(X = map_input, check_duplicates = FALSE, max_iter = 2000, normalize = FALSE,
                                  stop_lying_iter = 700, mom_switch_iter = 700,
                                  eta = round(nrow(map_input)/12),
                                  num_threads = ceiling(detectCores()/2))
    map <- map_calculate[["Y"]]
    colnames(map) <- c("tSNE1","tSNE2")
  }
  if(tolower(algorithm)=="umap") {
    if(tolower(language)=="r") {
      fcs_join_obj[["umap"]] <- list(coordinates = map,
                                          settings = list(language = "R", n_neighbors = 30, init = "spca", min_dist = 0.1,
                                                          n_threads = ceiling(detectCores()/2), verbose = TRUE))
    } else if(tolower(language)=="python") {
      fcs_join_obj[["umap"]] <- list(coordinates = map,
                                          settings = list(language = "Python", n_neighbors = 30, init = 'spectral',
                                                          min_dist = 0.1, low_memory = 'True', random_state = 123,
                                                          transform_seed = 123, verbose = 'True'))
    }
  } else if(tolower(algorithm)=="tsne") {
    if(tolower(language)=="r") {
      fcs_join_obj[["tsne"]] <- list(coordinates = map,
                                          settings = list(language = "R", check_duplicates = FALSE, max_iter = 2000,
                                                          normalize = FALSE, stop_lying_iter = 700, mom_switch_iter = 700,
                                                          eta = round(nrow(map_input)/12), num_threads = ceiling(detectCores()/2)))
    } else if(tolower(language)=="python") {
      stop("not supported (yet)")
    }
  }
  return(fcs_join_obj)
}
