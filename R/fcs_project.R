#' @title Project Marker Expression onto a 2D Embedding
#'
#' @description
#'   Creates a multi‐panel PDF of channel (parameter) expression projected onto
#'   a 2D reduction embedding (UMAP or tSNE). For each specified parameter,
#'   events are optionally trimmed of outliers, subsampled, and plotted as a
#'   rasterized scatter colored by expression intensity. Panels are arranged
#'   in 2×2 grids and written to disk.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join(), optionally augmented by
#'   fcs_batch_correction() and fcs_reduce_dimensions(), containing:
#'   - `data`: raw or transformed expression matrix (events × channels)
#'   - `batch_correction$data` (if present)
#'   - a `umap` or `tsne` element with `$coordinates` (events × 2)
#'
#' @param override_correction
#'   Logical; if `TRUE` (default), use `fcs_join_obj$data` even when
#'   batch correction exists. If `FALSE`, uses
#'   `fcs_join_obj$batch_correction$data`.
#'
#' @param reduction
#'   Character; embedding to use. Must be either `"UMAP"` (default) or `"tSNE"`.
#'
#' @param parameters
#'   Character vector of channel names to project. If `"all"` (default),
#'   all columns in the selected data matrix are plotted.
#'
#' @param outdir
#'   Character; path to an existing directory where the PDF will be saved.
#'   Defaults to `getwd()`.
#'
#' @param sample_size
#'   Integer; maximum number of events to sample per parameter before plotting
#'   (default 50000).
#'
#' @param point_size
#'   Numeric; point size for the rasterized scatter (`geom_point_rast`)
#'   (default 0.8).
#'
#' @param trim_outliers
#'   Logical; if `TRUE` (default), trims the lower and upper
#'   `trim_quantile` fraction of observations before sampling.
#'
#' @param trim_quantile
#'   Numeric; quantile threshold for outlier trimming (default 0.01).
#'
#' @param force_xlim
#'   Numeric vector of length 2 to fix the x‐axis limits, or `FALSE`
#'   (default) to use the data range.
#'
#' @param force_ylim
#'   Numeric vector of length 2 to fix the y‐axis limits, or `FALSE`
#'   (default) to use the data range.
#'
#' @details
#'   1. Determines which expression matrix to use based on
#'      `override_correction`.  
#'   2. Extracts UMAP or tSNE coordinates from
#'      `fcs_join_obj[[ tolower(reduction) ]][["coordinates"]]`.  
#'   3. Selects the requested `parameters` (all by default).  
#'   4. For each parameter:  
#'      - Trims outliers at the specified quantiles.  
#'      - Subsamples up to `sample_size` points.  
#'      - Creates a rasterized scatter plot colored by expression  
#'        (`scale_color_viridis(option="D")`).  
#'      - Honors any fixed `force_xlim`/`force_ylim`.  
#'   5. Arranges plots into pages of 2×2 panels via ggpubr::ggarrange().  
#'   6. Saves a timestamped PDF named  
#'      `<reduction>_parameter_projections_<YYYY-MM-DD_HHMMSS>.pdf`.  
#'
#' @return
#'   Invisibly returns `NULL`. The side effect is a multi‐page PDF written
#'   to `outdir`.
#'
#' @examples
#' \dontrun{
#'   # Assume joined, reduced and batch‐corrected object
#'   joined <- FCSimple::fcs_join(files)
#'   reduced <- FCSimple::fcs_reduce_dimensions(joined, method = "UMAP")
#'   corrected <- FCSimple::fcs_batch_correction(reduced)
#'
#'   # Project all channels onto UMAP
#'   FCSimple::fcs_project_parameters(
#'     corrected,
#'     reduction = "UMAP",
#'     outdir    = "~/results"
#'   )
#'
#'   # Project only CD3 and CD19
#'   FCSimple::fcs_project_parameters(
#'     corrected,
#'     parameters = c("CD3","CD19"),
#'     sample_size = 30000
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_reduce_dimensions, FCSimple::fcs_plot_reduction,
#'   ggplot2::ggplot, viridis::scale_color_viridis,
#'   ggrastr::geom_point_rast, ggpubr::ggarrange
#'
#' @importFrom ggplot2 ggplot aes theme_minimal theme xlim ylim ggtitle
#' @importFrom viridis scale_color_viridis
#' @importFrom ggrastr geom_point_rast
#' @importFrom ggpubr ggarrange
#' @importFrom gridExtra marrangeGrob
#' @export
fcs_project_parameters <- function(fcs_join_obj,
                                   override_correction = TRUE, 
                                   reduction = c("UMAP","tSNE"),
                                   parameters = "all",
                                   outdir = getwd(),
                                   sample_size = 50000,
                                   point_size = 0.8,
                                   point_alpha = 0.1, 
                                   trim_outliers = TRUE,
                                   trim_quantile = 0.01, 
                                   force_xlim = FALSE, 
                                   force_ylim = FALSE)
{
  require(ggplot2)
  require(viridis)
  require(ggpubr)
  require(ggrastr)

  if(length(reduction)!=1) {
    stop("error in argument 'reduction': use either 'UMAP' or 'tSNE'")
  }
  if("batch_correction" %in% names(fcs_join_obj)) {
    if(override_correction) {
      print("batch_correction found in fcs_join_obj and override_correction set to TRUE. Using fcs_join_obj[['data']] for projections. To use batch-corrected features, set 'override_correction' to FALSE.")
      join_data <- fcs_join_obj[["data"]]
    } else {
      print("batch_correction found in fcs_join_obj and override_correction set to FALSE. Using fcs_join_obj[['batch_correction']][['data']] for projections. To use original features, set 'override_correction' to TRUE.")
      join_data <- fcs_join_obj[["batch_correction"]][["data"]]
    }
  } else {
    join_data <- fcs_join_obj[["data"]]
  }
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
  intens_pl <- function(arg1, method = "color", pts = point_size, palpha = point_alpha, 
                        tr_out = trim_outliers, tr_q = trim_quantile) { # allow for later output as pch = 21 or similar with fill
    arg1 <- as.data.frame(arg1)
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
        plt <- ggplot(data = arg1, mapping = aes(x = UMAP1, y = UMAP2, color = col1))
      } else if(tolower(reduction)=="tsne") {
        plt <- ggplot(data = arg1, mapping = aes(x = tSNE1, y = tSNE2, color = col1))
      }
      plt <- plt +
        geom_point_rast(size = pts, pch = 19, alpha = palpha) +
        scale_color_viridis(option = "D", 
                            guide = guide_colorbar(title.position = 'top', 
                                                   frame.colour = "black",
                                                   ticks.colour = "black",
                                                   frame.linewidth = 0.4, draw.ulim = TRUE, draw.llim = TRUE,
                                                   label.theme = element_text(angle = 90, vjust = 0.5, size = 16), 
                                                   position = 'bottom', barwidth = 10, barheight = 1)) +
        ggtitle(cname1) +
        theme_bw(base_size = 18) + 
        theme(axis.title = element_blank(),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "bottom", 
              legend.spacing = unit(0, "pt"),
              legend.margin = margin(t = -5, r = 0, b = 20, l = 0),
              plot.margin = margin(t = 5, r = 5, b = 0, l = 5))
    } else if(method=="fill") {
      if(tolower(reduction)=="umap") {
        plt <- ggplot(data = arg1, mapping = aes(x = UMAP1, y = UMAP2, fill = col1))
      } else if(tolower(reduction)=="tsne") {
        plt <- ggplot(data = arg1, mapping = aes(x = tSNE1, y = tSNE2, fill = col1))
      }
      plt <- plt +
        geom_point_rast(size = pts, pch = 21, stroke = 0.05, alpha = palpha) +
        scale_fill_viridis(option = "D", 
                           guide = guide_colorbar(title.position = 'top', 
                                                  frame.colour = "black",
                                                  ticks.colour = "black",
                                                  frame.linewidth = 0.4, draw.ulim = TRUE, draw.llim = TRUE,
                                                  label.theme = element_text(angle = 90, vjust = 0.5, size = 16), 
                                                  position = 'bottom', barwidth = 8, barheight = 1)) +
        ggtitle(cname1) +
        theme_bw(base_size = 18) + 
        theme(axis.title = element_blank(),
              axis.text = element_blank(), 
              axis.ticks = element_blank(), 
              legend.title = element_blank(),
              plot.title = element_text(hjust = 0.5),
              legend.position = "bottom", 
              legend.spacing = unit(0, "pt"),
              legend.margin = margin(t = -5, r = 0, b = 20, l = 0),
              plot.margin = margin(t = 5, r = 5, b = 0, l = 5))
    }
    if(class(force_xlim)=='numeric') {
      plt <- plt + xlim(force_xlim)
    }
    if(class(force_ylim)=='numeric') {
      plt <- plt + ylim(force_ylim)
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
