fcs_plot_reduction <- function(fcs_join_obj, algorithm, reduction, point_alpha = 0.1, point_size = 1, outdir = getwd(),
                               split_factor = NA, internal_call = FALSE, anno_indices = NULL, keep_indices = NA,
                               pdf_dim = 10, png_dim = 1000, plotting_device = "pdf", annotate_text_size = 5,
                               title_size = 14, return_plot = TRUE, randomize_colors = FALSE, color_random_seed = 123,
                               color_clusters = TRUE, force_title = FALSE, sample_equally = TRUE)
{
  # use annotate_text_size = NA to produce a plot without cluster annotations
  require(ggplot2)
  require(shadowtext)
  require(ggrastr)

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  cluster_numbers <- as.numeric(as.character(fcs_join_obj[[tolower(algorithm)]][["clusters"]]))
  uclus <- unique(cluster_numbers)[order(unique(cluster_numbers))]

  # https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  colorby <- gg_color_hue(n = length(uclus))
  if(randomize_colors) {
    set.seed(color_random_seed)
    names(colorby) <- sample(x = uclus, size = length(uclus), replace = FALSE)
  } else {
    names(colorby) <- uclus
  }

  xval <- rep(NA,times=length(uclus)); names(xval) <- uclus; yval <- xval
  for(i in 1:length(xval)) {
    xval[i] <- median(reduction_coords[,1][which(cluster_numbers==as.numeric(names(xval)[i]))])
    yval[i] <- median(reduction_coords[,2][which(cluster_numbers==as.numeric(names(yval)[i]))])
  }
  plt_input <- cbind(reduction_coords,data.frame(cluster = cluster_numbers))
  plt_input$cluster <- factor(plt_input$cluster)
  if(!internal_call) {
    pl_fun <- function(plin, ptalpha = point_alpha, xanno = xval,
                       yanno = yval, sizeanno = annotate_text_size, ftitle = force_title,
                       color_clus = color_clusters, ptsize = point_size)
    {
      cnamex <- colnames(plin)[1]; cnamey <- colnames(plin)[2]
      colnames(plin)[1:2] <- c("valx","valy") # bypass aes_string, which is now deprecated
      if(color_clus) {
        mypl <- ggplot(data = plin, mapping = aes(x = valx,
                                                  y = valy,
                                                  color = cluster)) +
          ggrastr::geom_point_rast(alpha = ptalpha, size = ptsize) +
          scale_color_manual(values = colorby) +
          labs(x = cnamex, y = cnamey)
      } else {
        mypl <- ggplot(data = plin, mapping = aes(x = valx,
                                                  y = valy)) +
          ggrastr::geom_point_rast(alpha = ptalpha, size = ptsize) +
          labs(x = cnamex, y = cnamey)
      }
      if(!is.na(sizeanno)) {
        mypl <- mypl +
          annotate("shadowtext", x = xanno, y = yanno, label = names(xanno), size = sizeanno)
      }
      mypl <- mypl + theme_void() +
        theme(legend.position = "none")
      if(ftitle) {
        mypl <- mypl + ggtitle(colnames(plin)[ncol(plin)]) + theme(plot.title = element_text(hjust = 0.5, size = title_size))
      }
      return(mypl)
    }
    if(is.na(split_factor[1])) {
      plt_reduction <- pl_fun(plin = plt_input)
    } else {
      # split the plot by split_factor where split_factor is a vector of length == nrow(reduction) that gives identity to the reduction rows
      split_reduction <- split(x = plt_input, f = factor(split_factor))
      for(i in 1:length(split_reduction)) {
        split_reduction[[i]] <- cbind(as.data.frame(split_reduction[[i]]), data.frame(var1 = rep(1,nrow(split_reduction[[i]]))))
        colnames(split_reduction[[i]])[ncol(split_reduction[[i]])] <- names(split_reduction)[i]
      }
      if(sample_equally) {
        sample_size <- min(as.numeric(sapply(split_reduction,nrow)))
        for(i in 1:length(split_reduction)) {
          if(nrow(split_reduction[[i]])>sample_size) {
            set.seed(123)
            split_reduction[[i]] <- split_reduction[[i]][sample(1:nrow(split_reduction[[i]]),sample_size,replace=F),]
          }
        }
      }
      outplots <- lapply(X = split_reduction, pl_fun, ftitle = force_title)
      plotnrow <- ifelse(length(outplots)<=4,1,floor(sqrt(length(outplots))))
      plotncol <- ifelse(length(outplots)<=4,length(outplots),ceiling(sqrt(length(outplots))))
      plt_reduction <- ggpubr::ggarrange(plotlist = outplots, nrow = plotnrow, ncol = plotncol)
    }
    fname <- paste0(outdir,"/",tolower(algorithm),"_",tolower(reduction),"_labeled_",
                    strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
  } else {
    plt_reduction <- ggplot(data = plt_input[-keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                       y = colnames(reduction_coords)[2])) +
      ggrastr::geom_point_rast(alpha = point_alpha, color = "grey") +
      ggrastr::geom_point_rast(data = plt_input[keep_indices,], mapping = aes_string(x = colnames(reduction_coords)[1],
                                                                                        y = colnames(reduction_coords)[2]),
                               alpha = point_alpha, color = "red") +
      annotate("shadowtext", x = xval, y = yval, label = names(xval), size = annotate_text_size) +
      theme_void() +
      theme(legend.position = "none")
    fname <- paste0(outdir,"/islands_selected_for_by_dbscan_",
                    strftime(Sys.time(),"%Y-%m-%d_%H%M%S"))
  }
  if(return_plot) {
    return(plt_reduction)
  } else {
    if(plotting_device=="pdf") {
      ggsave(filename = paste0(fname,".pdf"),
             plot = plt_reduction, device = "pdf", width = pdf_dim, height = pdf_dim,
             units = "in", dpi = 900)
    } else if(plotting_device=="png"){
      ggsave(filename = paste0(fname,".png"),
             plot = plt_reduction, device = "png", width = png_dim, height = png_dim,
             units = "in", dpi = 900)
    }
  }
}
