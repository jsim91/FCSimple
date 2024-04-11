fcs_plot_reduction_density <- function(fcs_join_obj, reduction = "UMAP", n_kde = 200,
                                       xlimits = "auto", ylimits = "auto",
                                       z_inflation = 200, text_size_factor = 1,
                                       title_string = "my title", n_color_steps = 15,
                                       RColorBrewer_pal = "Spectral",
                                       plot_border_thickness = 1)
{
  require(ggplot2)
  require(MASS)
  require(reshape2)
  require(scales)
  require(RColorBrewer)
  require(metR)

  if(xlimits[1]=="auto") {
    xlimits = range(fcs_join_obj$umap$coordinates[,1])
  }
  if(ylimits[1]=="auto") {
    ylimits = range(fcs_join_obj$umap$coordinates[,2])
  }

  if(length(reduction)!=1) {
    stop("error in argument 'reduction': Use either tSNE or UMAP. Reduction must be present in 'fcs_join_obj'.")
  }

  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  reduction_names <- colnames(reduction_coords)[1:2]
  colnames(reduction_coords)[1:2] <- c("dimx","dimy")
  dim1_range <- xlimits; dim2_range <- ylimits

  if(length(cluster_algorithm)!=1) {
    cluster_algorithm <- cluster_algorithm[1]
  }
  if(annotate_clusters) {
    obj_clusters <- fcs_join_obj[[tolower(cluster_algorithm)]][["clusters"]]
    unique_clus <- unique(obj_clusters); unique_clus <- unique_clus[order(as.numeric(unique_clus))]
    clusx <- rep(NA,length(unique_clus)); names(clusx) <- unique_clus; clusy <- clusx
    for(i in 1:length(clusx)) {
      clusx[i] <- median(reduction_coords[,1][which(obj_clusters==as.numeric(names(clusx)[i]))])
      clusy[i] <- median(reduction_coords[,2][which(obj_clusters==as.numeric(names(clusx)[i]))])
    }
  }

  grp1_kde2d <- MASS::kde2d(x = reduction_coords[,1], y = reduction_coords[,2], lims = c(dim1_range, dim2_range), n = n_kde)
  diff12 = grp1_kde2d
  diff12$z <- grp1_kde2d$z
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y
  diff12.m = melt(diff12$z, id.var=rownames(diff12))
  colnames(diff12.m)[3] <- "z"

  craftbrewer_pal <- function (type = "seq", palette = 1, direction = 1)
  {
    pal <- scales:::pal_name(palette, type)
    force(direction)
    function(n) {
      n_max_palette <- RColorBrewer:::maxcolors[names(RColorBrewer:::maxcolors) == palette]

      if (n < 3) {
        pal <- suppressWarnings(RColorBrewer::brewer.pal(n, pal))
      } else if (n > n_max_palette){
        rlang::warn(paste(n, "colours used, but", palette, "has only",
                          n_max_palette, "- New palette created based on all colors of",
                          palette))
        n_palette <- RColorBrewer::brewer.pal(n_max_palette, palette)
        colfunc <- grDevices::colorRampPalette(n_palette)
        pal <- colfunc(n)
      }
      else {
        pal <- RColorBrewer::brewer.pal(n, pal)
      }
      pal <- pal[seq_len(n)]
      if (direction == -1) {
        pal <- rev(pal)
      }
      # Force the first color to be white
      pal[1] <- "white"
      pal
    }
  }

  scale_fill_craftfermenter <- function(..., type = "seq", palette = 1, direction = -1, na.value = "white", guide = "coloursteps", aesthetics = "fill") {
    type <- match.arg(type, c("seq", "div", "qual"))
    if (type == "qual") {
      warn("Using a discrete colour palette in a binned scale.\n  Consider using type = \"seq\" or type = \"div\" instead")
    }
    binned_scale(aesthetics, "fermenter", ggplot2:::binned_pal(craftbrewer_pal(type, palette, direction)), na.value = na.value, guide = guide, ...)
  }

  diff12.m.copy <- diff12.m
  mybreaks <- seq(from = min(200*diff12.m$z), to = max(200*diff12.m$z), length.out = n_color_steps)

  plt <- ggplot() +
    metR::geom_contour_fill(data = diff12.m, mapping = aes(Var1, Var2, z = 200*z)) +
    scale_fill_craftfermenter(
      breaks = mybreaks,
      palette = RColorBrewer_pal,
      guide = guide_colorsteps(
        frame.colour = "black",
        ticks.colour = "black",
        barwidth=20,
        )
    ) +
    geom_contour(data = diff12.m.copy, mapping = aes(Var1, Var2, z = 200*z), color = "black") +
    theme(legend.position = "bottom") +
    coord_cartesian(xlim = dim1_range, ylim = dim2_range, expand = FALSE) +
    labs(x = reduction_names[1], y = reduction_names[2], title = title_string)
  if(!is.na(plot_border_thickness[1])) {
    plt <- plt + theme_bw() +
      theme(axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 28*text_size_factor),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.title = element_blank(),
            legend.position = "none",
            panel.border = element_rect(color = "black", size = plot_border_thickness))
  } else {
    plt <- plt + theme(axis.text = element_blank(),
                       plot.title = element_text(hjust = 0.5, size = 28*text_size_factor),
                       axis.title = element_blank(),
                       axis.ticks = element_blank(),
                       legend.title = element_blank(),
                       legend.position = "none")
  }
  return(plt)
}
