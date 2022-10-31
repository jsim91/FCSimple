fcs_plot_distribution <- function(fcs_join_obj,
                                  separate_by = c("none", "date", "cluster"),
                                  plot_algorithm = c("leiden","flowsom","louvain","phenograph"),
                                  outdir = getwd(),
                                  plot_palette = NULL)
{
  require(ggpubr)
  require(ggplot2)
  require(ggridges)

  obj_data <- fcs_join_obj[["data"]]
  if(separate_by=="date") {
    if(!"run_date" %in% names(fcs_join_obj)){
      print("Unable to find run date. Using separate_by = 'none' instead.")
    } else {
      batch <- fcs_join_obj[["run_date"]]
      obj_data$date <- batch
    }
  }
  if("date" %in% colnames(obj_data)) {
    data_split <- vector("list", length = ncol(obj_data)-1)
    names(data_split) <- colnames(obj_data)[1:(ncol(obj_data)-1)]
    for(i in 1:length(data_split)) {
      data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])],
                                    batch = obj_data[,"batch"])
    }
  } else {
    data_split <- vector("list", length = ncol(obj_data))
    names(data_split) <- colnames(obj_data)
    for(i in 1:length(data_split)) {
      data_split[[i]] <- data.frame(val1 = obj_data[,which(colnames(obj_data)==names(data_split)[i])])
      colnames(data_split[[i]]) <- colnames(obj_data)[i]
    }
  }
  if(tolower(separate_by) == "cluster") {
    if(!tolower(plot_algorithm) %in% c("leiden","flowsom","louvain","phenograph")) {
      stop(paste0("error in argument 'plot_algorithm': clusters not found for ",tolower(plot_algorithm)," not found in 'fcs_join_obj'"))
    }
    for(i in 1:length(data_split)) {
      data_split[[i]]$cluster <- factor(fcs_join_obj[[tolower(plot_algorithm)]]$clusters)
    }
  }
  if(separate_by == "none") {
    plot_set <- lapply(X = data_split, FUN = plot_none)
    ggsave(filename = paste0(outdir,"/panel_distributions_concatenation_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
           plot = ggpubr::ggarrange(plotlist = plot_set, nrow = ceiling(sqrt(length(plot_set))),
                                    ncol = ceiling(sqrt(length(plot_set)))), device = "pdf",
           width = ceiling(sqrt(length(plot_set)))*2.5, height = ceiling(sqrt(length(plot_set)))*2.5, units = "in",
           dpi = 900, limitsize = FALSE)
  } else if(separate_by == "date") {
    plot_set <- lapply(X = data_split, FUN = plot_date)
    ggsave(filename = paste0(outdir,"/panel_distributions_by_batch_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
           plot = ggpubr::ggarrange(plotlist = plot_set, nrow = ceiling(sqrt(length(plot_set))),
                                    ncol = ceiling(sqrt(length(plot_set))), legend = "bottom", common.legend = TRUE),
           device = "pdf", width = ceiling(sqrt(length(plot_set)))*2.5, height = ceiling(sqrt(length(plot_set)))*2.5,
           units = "in", dpi = 900, limitsize = FALSE)
  } else if(separate_by == "cluster") {
    plot_set <- lapply(X = data_split, FUN = plot_cluster)
    ggsave(filename = paste0(outdir,"/panel_distributions_by_cluster_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
           plot = ggpubr::ggarrange(plotlist = plot_set, nrow = ceiling(sqrt(length(plot_set))),
                                    ncol = ceiling(sqrt(length(plot_set))), legend = "bottom", common.legend = TRUE),
           device = "pdf", width = ceiling(sqrt(length(plot_set)))*3, height = ceiling(sqrt(length(plot_set)))*9,
           units = "in", dpi = 900, limitsize = FALSE)
  }
}

plot_none <- function(input_data)
{
  xval <- density(input_data[,1])$x; yval <- density(input_data[,])$y
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

plot_date <- function(input_data)
{
  capture_channel <- colnames(input_data)[1]
  colnames(input_data)[1] <- "tmp"
  plt1 <- ggplot(data = input_data, mapping = aes(x = tmp, group = batch)) +
    geom_density(lwd = 0.3) +
    theme_minimal() +
    ggtitle(capture_channel) +
    theme(axis.text.y = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5))
  return(plt1)
}

plot_cluster <- function(input_data)
{
  capture_channel <- colnames(input_data)[1]
  colnames(input_data)[1] <- "tmp"
  plt1 <- ggplot(data = input_data, mapping = aes(x = tmp, y = cluster)) +
    geom_density_ridges(lwd = 0.3) +
    theme_minimal() +
    ggtitle(capture_channel) +
    theme(axis.text.y = element_text(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  return(plt1)
}

