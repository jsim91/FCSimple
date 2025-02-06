fcs_plot_distribution <- function(fcs_join_obj,
                                  override_correction = TRUE, 
                                  separate_by = c("none", "date", "cluster"),
                                  plot_algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                                  outdir = getwd(),
                                  plot_palette = NULL,
                                  rm_zero = FALSE,
                                  trim_quantile = NULL,
                                  add_timestamp = TRUE)
{
  require(ggpubr)
  require(ggplot2)
  require(ggridges)

  use_rep <- tolower(use_rep)
  if('batch_correction' %in% names(fcs_join_obj)) {
    if("override_correction") {
      print("Batch_correction found in fcs_join_obj list with 'override_correction' set to TRUE. Using fcs_join_obj[['data']] for plotting. To use batch-corrected features, set 'override_correction' to FALSE.")
      obj_data <- fcs_join_obj[['data']]
    } else {
      print("Batch_correction found in fcs_join_obj list with 'override_correction' set to FALSE. Using fcs_join_obj[['batch_correction']][['data']] for plotting. To use uncorrected features, set 'override_correction' to TRUE")
      obj_data <- fcs_join_obj[['batch_correction']][['data']]
    }
  } else {
    print("batch_correction not found in fcs_join_obj list. Proceeding with fcs_join_obj[['data']].")
    obj_data <- as.data.frame(fcs_join_obj[["data"]])
  }
  if(!'collection_instrument' %in% names(fcs_join_obj)) {
    stop("'fcs_join_obj[['collection_instrument']] not found; do fcs_update() on your object and specify what instrument_type the data was collected with first, either 'cytof' for mass cytometry or 'flow' for fluorescence cytometry")
  }
  if(tolower(separate_by)=="date") {
    if(!"run_date" %in% names(fcs_join_obj)){
      print("Unable to find run date. Using separate_by = 'none' instead.")
    } else {
      batch <- fcs_join_obj[["run_date"]]
      obj_data$date <- batch
    }
    data_split <- vector("list", length = ncol(obj_data)-1)
    names(data_split) <- colnames(obj_data)[1:(ncol(obj_data)-1)]
    for(i in 1:length(data_split)) {
      data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])],
                                    batch = obj_data[,"date"])
      colnames(data_split[[i]])[1] <- colnames(obj_data)[i]
    }
  } else if(tolower(separate_by) == "cluster") {
    data_split <- vector("list", length = ncol(obj_data))
    names(data_split) <- colnames(obj_data)
    for(i in 1:length(data_split)) {
      data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])])
      colnames(data_split[[i]])[1] <- colnames(obj_data)[i]
    }
    if(!tolower(plot_algorithm) %in% c("leiden","flowsom","louvain","phenograph","git")) {
      stop(paste0("error in argument 'plot_algorithm': clusters not found for ",tolower(plot_algorithm)," not found in 'fcs_join_obj'"))
    }
    for(i in 1:length(data_split)) {
      data_split[[i]]$cluster <- factor(fcs_join_obj[[tolower(plot_algorithm)]]$clusters)
    }
  }
  if(tolower(separate_by) == "none") {
    data_split <- vector("list", length = ncol(obj_data))
    names(data_split) <- colnames(obj_data)
    for(i in 1:length(data_split)) {
      data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])])
      colnames(data_split[[i]])[1] <- colnames(obj_data)[i]
    }
    plot_set <- lapply(X = data_split, FUN = plot_none, rm0 = rm_zero, trq = trim_quantile, instr_type = fcs_join_obj[['collection_instrument']])
    if(add_timestamp) {
      fname <- paste0(outdir,"/panel_distributions_concatenation_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf")
    } else {
      fname <- paste0(outdir,"/panel_distributions_concatenation.pdf")
    }
    ggsave(filename = fname,
           plot = ggpubr::ggarrange(plotlist = plot_set, nrow = ceiling(sqrt(length(plot_set))),
                                    ncol = ceiling(sqrt(length(plot_set)))), device = "pdf",
           width = ceiling(sqrt(length(plot_set)))*2.5, height = ceiling(sqrt(length(plot_set)))*2.5, units = "in",
           dpi = 900, limitsize = FALSE)
  } else if(tolower(separate_by) == "date") {
    if(!"batch" %in% colnames(data_split[[1]])) {
      stop("error in argument 'separate_by': no run date found, cannot plot by run date/batch")
    }
    plot_set <- lapply(X = data_split, FUN = plot_date, rm0 = rm_zero, trq = trim_quantile, instr_type = fcs_join_obj[['collection_instrument']])
    if(add_timestamp) {
      fname <- paste0(outdir,"/panel_distributions_by_batch_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf")
    } else {
      fname <- paste0(outdir,"/panel_distributions_by_batch.pdf")
    }
    ggsave(filename = fname,
           plot = ggpubr::ggarrange(plotlist = plot_set, nrow = ceiling(sqrt(length(plot_set))),
                                    ncol = ceiling(sqrt(length(plot_set))), legend = "bottom", common.legend = TRUE),
           device = "pdf", width = ceiling(sqrt(length(plot_set)))*2.5, height = ceiling(sqrt(length(plot_set)))*2.5,
           units = "in", dpi = 900, limitsize = FALSE)
  } else if(tolower(separate_by) == "cluster") {
    plot_set <- lapply(X = data_split, FUN = plot_cluster, rm0 = rm_zero, trq = trim_quantile, instr_type = fcs_join_obj[['collection_instrument']])
    if(add_timestamp) {
      fname <- paste0(outdir,"/panel_distributions_by_cluster_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf")
    } else {
      fname <- paste0(outdir,"/panel_distributions_by_cluster.pdf")
    }
    ggsave(filename = fname,
           plot = ggpubr::ggarrange(plotlist = plot_set, nrow = ceiling(sqrt(length(plot_set))),
                                    ncol = ceiling(sqrt(length(plot_set))), legend = "bottom", common.legend = TRUE),
           device = "pdf", width = ceiling(sqrt(length(plot_set)))*3, height = ceiling(sqrt(length(plot_set)))*9,
           units = "in", dpi = 900, limitsize = FALSE)
  }
}

plot_none <- function(input_data, rm0, trq, instr_type)
{
  cache_colname <- colnames(input_data)[1]
  if(instr_type=="flow") {
    if(rm0) {
      warning("remove 0s set to TRUE while instrument_type is 'flow'; leaving 0 values.")
    }
    densdat <- input_data[,1]
    if(!is.null(trq)) {
      if(length(trq)!=2) {
        stop("'trim_quantile' should be a vector of length 2")
      } 
      trim_q <- as.numeric(quantile(x = densdat, probs = trq))
      drop_index <- union(which(densdat>max(trim_q)),which(densdat<min(trim_q)))
      densdat <- densdat[-drop_index]
    }
    xval <- density(densdat)$x; yval <- density(densdat)$y
  }
  if(instr_type=="cytof") {
    if(rm0) {
      rm0_ind <- which(input_data[,1]==0)
      if(length(rm0_ind)!=0) {
        densdat <- input_data[-which(input_data[,1]==0),1]
      } else {
        densdat <- input_data[,1]
      }
      if(!is.null(trq)) {
        if(length(trq)!=1) {
          stop("'trim_quantile' should be a singlet numeric value")
        } 
        trim_q <- as.numeric(quantile(x = densdat, probs = trq))
        densdat <- densdat[-which(densdat>trim_q)]
      }
      xval <- density(densdat)$x; yval <- density(densdat)$y
    } else {
      if(!is.null(trq)) {
        trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
        densdat <- input_data[-which(input_data[,1]>trim_q),]
      } else {
        densdat <- input_data[,1]
      }
      xval <- density(densdat)$x; yval <- density(densdat)$y
    }
  }
  dens_data <- data.frame(xval = xval, yval = yval)
  plt1 <- ggplot(data = dens_data, mapping = aes(x = xval, y = yval)) +
    geom_ribbon(mapping = aes(ymin = 0, ymax = yval), color = "grey", alpha = 0.5) +
    geom_line(lwd = 0.3) +
    theme_minimal() +
    ggtitle(colnames(input_data)[1]) +
    theme(axis.text.y = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  return(plt1)
}

plot_date <- function(input_data, rm0, trq, instr_type)
{
  capture_channel <- colnames(input_data)[1]
  colnames(input_data)[1] <- "tmp"
  if(instr_type=="flow") {
    if(rm0) {
      warning("remove 0s set to TRUE while instrument_type is 'flow'; leaving 0 values.")
    }
    if(!is.null(trq)) {
      if(length(trq)!=2) {
          stop("'trim_quantile' should be a vector of length 2")
        }
      trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
      drop_index <- union(which(input_data[,1]>max(trim_q)),which(input_data[,1]<min(trim_q)))
      input_data <- input_data[-drop_index,]
    }
  }
  if(instr_type=="cytof") {
    if(rm0) {
      rm0_ind <- which(input_data[,1]==0)
      if(length(rm0_ind)!=0) {
        input_data <- input_data[-which(input_data[,1]==0),]
      }
      if(!is.null(trq)) {
        if(length(trq)!=1) {
          stop("'trim_quantile' should be a singlet numeric value")
        } 
        trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
        input_data <- input_data[-which(input_data[,1]>trim_q),]
      }
    } else {
      if(!is.null(trq)) {
        if(length(trq)!=1) {
          stop("'trim_quantile' should be a singlet numeric value")
        } 
        trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
        input_data <- input_data[-which(input_data[,1]>trim_q),]
      }
    }
  }
  plt1 <- ggplot(data = input_data, mapping = aes(x = tmp, group = batch,
                                                  color = batch, fill = batch)) +
    geom_density(size = 0.75, alpha = 0.1) +
    theme_minimal() +
    ggtitle(capture_channel) +
    theme(axis.text.y = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  return(plt1)
}

plot_cluster <- function(input_data, rm0, trq, instr_type)
{
  capture_channel <- colnames(input_data)[1]
  colnames(input_data)[1] <- "tmp"
  if(instr_type=="flow") {
    if(rm0) {
      warning("remove 0s set to TRUE while instrument_type is 'flow'; leaving 0 values.")
    }
    if(!is.null(trq)) {
      if(length(trq)!=2) {
          stop("'trim_quantile' should be a vector of length 2")
        }
      trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
      drop_index <- union(which(input_data[,1]>max(trim_q)),which(input_data[,1]<min(trim_q)))
      input_data <- input_data[-drop_index,]
    }
  }
  if(instr_type=="cytof") {
    if(rm0) {
      rm0_ind <- which(input_data[,1]==0)
      if(length(rm0_ind)!=0) {
        input_data <- input_data[-which(input_data[,1]==0),]
      }
      if(!is.null(trq)) {
        trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
        input_data <- input_data[-which(input_data[,1]>trim_q),]
      }
    } else {
      if(!is.null(trq)) {
        trim_q <- as.numeric(quantile(x = input_data[,1], probs = trq))
        input_data <- input_data[-which(input_data[,1]>trim_q),]
      }
    }
  }
  plt1 <- ggplot(data = input_data, mapping = aes(x = tmp, y = cluster)) +
    geom_density_ridges(lwd = 0.3) +
    theme_minimal() +
    ggtitle(capture_channel) +
    theme(axis.text.y = element_text(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  return(plt1)
}
