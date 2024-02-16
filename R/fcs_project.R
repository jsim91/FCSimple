fcs_project_parameters <- function(fcs_join_obj,
                                   reduction = c("UMAP","tSNE"),
                                   parameters = "all",
                                   outdir = getwd(),
                                   sample_size = 50000,
                                   point_size = 0.8,
                                   trim_outliers = TRUE,
                                   trim_quantile = 0.01)
{
  require(ggplot2)
  require(viridis)
  require(ggpubr)
  require(ggrastr)

  if(length(reduction)!=1) {
    stop("error in argument 'reduction': use either 'UMAP' or 'tSNE'")
  }
  join_data <- fcs_join_obj[["data"]]
  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  if(nrow(join_data)<1) {
    stop("error in 'reduction': unable to find specified reduction")
  }
  if(parameters[1] == "all") {
    include_params <- colnames(join_data)
  } else {
    target_params <- colnames(join_data)[which(colnames(join_data) %in% parameters)]
    if(length(target_params)==0) {
      stop("error in argument 'parameters': none of the parameters requested were found")
    } else {
      include_params <- target_params
    }
  }
  intens_list <- vector("list", length = length(include_params)); names(intens_list) <- include_params
  for(i in 1:length(include_params)) {
    tmp_param <- include_params[i]
    intens_list[[i]] <- cbind(join_data[,tmp_param],reduction_coords); colnames(intens_list[[i]])[1] <- names(intens_list)[i]
  }
  intens_pl <- function(arg1, method = "color", pts = point_size,
                        tr_out = trim_outliers, tr_q = trim_quantile) { # allow for later output as pch = 21 or similar with fill
    if(tr_out) {
      quant_val <- quantile(x = arg1[,1], probs = c(tr_q, 1-tr_q))
      arg1 <- arg1[intersect(which(arg1[,1]>min(quant_val)), which(arg1[,1]<max(quant_val))),]
    }
    if(nrow(arg1)>sample_size){
      arg1 <- arg1[sample(1:nrow(arg1),size=sample_size,replace=F),]
    }
    cname1 <- colnames(arg1)[1]
    colnames(arg1)[1] <- "col1"
    if(method=="color") {
      if(tolower(reduction)=="umap"){
        plt <- ggplot(data = arg1, mapping = aes_string(x = "UMAP1", y = "UMAP2", color = colnames(arg1)[1]))
      } else if(tolower(reduction)=="tsne") {
        plt <- ggplot(data = arg1, mapping = aes_string(x = "tSNE1", y = "tSNE2", color = colnames(arg1)[1]))
      }
      plt <- plt +
        geom_point_rast(size = pts, pch = 19, alpha = 0.1) +
        scale_color_viridis(option = "D") +
        ggtitle(cname1) +
        theme_minimal() +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "bottom")
    } else if(method=="fill") {
      if(tolower(reduction)=="umap") {
        plt <- ggplot(data = arg1, mapping = aes_string(x = "UMAP1", y = "UMAP2", fill = colnames(arg1)[1]))
      } else if(tolower(reduction)=="tsne") {
        plt <- ggplot(data = arg1, mapping = aes_string(x = "tSNE1", y = "tSNE2", fill = colnames(arg1)[1]))
      }
      plt <- plt +
        geom_point_rast(size = pts, pch = 21, stroke = 0.05) +
        scale_fill_viridis(option = "D") +
        ggtitle(cname1) +
        theme_minimal() +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "bottom")
    }
    return(plt)
  }
  intens_plots <- lapply(intens_list, intens_pl)

  void_plot <- ggplot(data = cars, aes(x=speed,y=dist)) +
    geom_point(alpha = 0) +
    theme_void()

  arrange_length <- ceiling(length(intens_list)/4)
  for(i in 1:(arrange_length*4)) {
    if(i>length(intens_plots)) {
      intens_plots[[i]] <- void_plot
    }
  }

  arranged_list <- vector("list", length = arrange_length)
  low_index <- seq(from = 1, to = arrange_length*4, by = 4)
  high_index <- seq(from = 4, to = arrange_length*4, by = 4)

  for(i in 1:length(arranged_list)) {
    arranged_list[[i]] <- ggpubr::ggarrange(plotlist = intens_plots[low_index[i]:high_index[i]],
                                            nrow = 2, ncol = 2)
  }

  ggsave(filename = paste0(ifelse(tolower(reduction)=="umap","UMAP","tSNE"),"_parameter_projections_",
                           strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
         plot = gridExtra::marrangeGrob(grobs = arranged_list, nrow=1, ncol=1, top = ""),
         device = "pdf", path = outdir, width = 12, height = 12, units = "in", dpi = 900)
}
