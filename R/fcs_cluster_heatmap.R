#' @title Generate Cluster Heatmap for Flow Cytometry Data
#'
#' @description
#'   Computes and visualizes median‐scaled expression of each cluster across
#'   parameters using ComplexHeatmap. Selects raw or batch‐corrected data,
#'   aggregates by cluster, and stores both the heatmap object and tile data.
#'
#' @param fcs_join_obj
#'   A list returned by FCSimple::fcs_join() and FCSimple::fcs_cluster(),
#'   containing at least:
#'   - `data`: numeric matrix of events × channels
#'   - `<algorithm>$clusters`: vector of cluster IDs per event
#'   - optionally `batch_correction$data` if correction was applied
#'
#' @param algorithm
#'   Character; name of clustering results in `fcs_join_obj` to use
#'   (e.g. `"leiden"`, `"flowsom"`, etc.).
#'
#' @param include_parameters
#'   Character vector of channel names to include (default `"all"`).
#'   If `"all"`, uses all columns of the expression matrix.
#'
#' @param include_clusters
#'   Character vector of cluster IDs to include (default `"all"`).
#'   If `"all"`, includes every cluster.
#'
#' @param heatmap_color_palette
#'   Character vector of colors (length ≥ 2) for the heatmap palette.
#'   Default uses a reversed “RdYlBu” from RColorBrewer.
#'
#' @param transpose_heatmap
#'   Logical; if `TRUE`, transpose the heatmap matrix before plotting.
#'   Default `FALSE`.
#'
#' @param cluster_row
#'   Logical; if `TRUE`, apply hierarchical clustering to rows.
#'   Default `TRUE`.
#'
#' @param cluster_col
#'   Logical; if `TRUE`, apply hierarchical clustering to columns.
#'   Default `TRUE`.
#'
#' @param override_correction
#'   Logical; if `TRUE`, always use raw data even if batch correction
#'   is present. Default `TRUE`.
#'
#' @param return_heatmap_data
#'   Logical; if `TRUE`, return the matrix of median‐scaled values
#'   invisibly instead of adding the heatmap to `fcs_join_obj`.
#'   Default `FALSE`.
#'
#' @param heatmap_linewidth
#'   Numeric; border line width for heatmap cells. Default `0.5`.
#'
#' @details
#' - Chooses raw or batch‐corrected data based on `override_correction` and
#'   presence of `fcs_join_obj$batch_correction$data`.
#' - Scales expression values to [0,1] per channel using
#'   CATALYST:::.scale_exprs.
#' - Computes median expression for each cluster × parameter.
#' - Builds a ComplexHeatmap object with cluster‐size annotations.
#' - Stores the result under
#'   `fcs_join_obj[[paste0(tolower(algorithm), "_heatmap")]]`:
#'     - `heatmap`: the Heatmap object
#'     - `heatmap_tile_data`: the numeric matrix used
#'     - `population_size`: cluster event counts
#'     - `rep_used`: “with batch correction” or “without batch correction”
#' - Appends a timestamped entry to `object_history`.
#'
#' @return
#'   Invisibly returns the updated `fcs_join_obj` with a new element
#'   `<algorithm>_heatmap` as described above. If
#'   `return_heatmap_data = TRUE`, returns only the heatmap matrix.
#'
#' @examples
#' \dontrun{
#'   joined <- FCSimple::fcs_join(list(ff1, ff2))
#'   clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")
#'
#'   # Generate and store heatmap object
#'   out <- FCSimple::fcs_cluster_heatmap(
#'     clustered,
#'     algorithm = "leiden",
#'     include_parameters = c("CD3","CD4","CD8"),
#'     override_correction = FALSE
#'   )
#'
#'   # Just get the tile matrix
#'   mat <- FCSimple::fcs_cluster_heatmap(
#'     clustered,
#'     algorithm = "leiden",
#'     return_heatmap_data = TRUE
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_cluster, FCSimple::fcs_plot_heatmap
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid unit gpar
#' @export
fcs_cluster_heatmap <- function(fcs_join_obj, algorithm, include_parameters = "all", include_clusters = "all",
                                heatmap_color_palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                transpose_heatmap = FALSE, cluster_row = TRUE, cluster_col = TRUE,
                                override_correction = TRUE, return_heatmap_data = FALSE,
                                heatmap_linewidth = 0.5, row_text_size = 13, column_text_size = 12,
                                legend_text_size = 11, scaling_quantile = 0.01)
{
  if(!tolower(algorithm) %in% names(fcs_join_obj)) {
    stop("error in argument 'algorithm': algorithm not found in fcs_join_obj. Try 'View(fcs_join_obj)'")
  }

  require(CATALYST)
  require(ComplexHeatmap)
  require(circlize)
  require(grid)

  if('batch_correction' %in% names(fcs_join_obj)) {
    if(override_correction) {
      cordat <- FALSE
      print("batch_correction found in fcs_join_obj list but using uncorrected data for heatmap because 'override_correction' was set to TRUE. Set 'override_correction' to FALSE to use corrected data.")
      heatmap_data <- fcs_join_obj[["data"]]
    } else {
      cordat <- TRUE
      print("batch_correction found in fcs_join_obj list. Using fcs_join_obj[['batch_correction']][['data']] for heatmap. Set 'override_correction' to TRUE to use uncorrected data.")
      heatmap_data <- fcs_join_obj[['batch_correction']][['data']]
      if(!'object_history' %in% names(fcs_join_obj)) {
        print("'object_history' not found in fcs_join_obj. Consider running FCSimple::fcs_update() on the object. Trying to use fcs_join_obj[['batch_correction']][['data']] for heatmap.")
      } else {
        if(sum(grepl(pattern = "pca", x = fcs_join_obj[['object_history']]))!=0) {
          print("It seems pca was run on this object. Using fcs_join_obj[['batch_correction']][['data']] will plot the median scaled expression scores of the corrected PCs.")
        }
      }
    }
  } else {
    cordat <- FALSE
    print("batch_correction not found in fcs_join_obj list. Using fcs_join_obj[['data']] for heatmap.")
    heatmap_data <- fcs_join_obj[["data"]]
  }
  event_source <- fcs_join_obj[["source"]]
  cluster_numbers <- as.character(fcs_join_obj[[tolower(algorithm)]][["clusters"]])
  if(include_parameters[1]=="all") {
    include_channels <- colnames(heatmap_data)
  } else {
    include_channels <- include_parameters
  }

  scaled.global <- CATALYST:::.scale_exprs(t(heatmap_data[,include_channels]), 1, q = scaling_quantile)
  global.t <- t(scaled.global)
  backend.matrix <- matrix(data=NA,nrow=length(unique(cluster_numbers)),ncol=ncol(global.t))
  row.names(backend.matrix) <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  colnames(backend.matrix) <- colnames(global.t)
  for(i in 1:nrow(backend.matrix)) {
    get.clus <- which(cluster_numbers==row.names(backend.matrix)[i])
    for(j in 1:ncol(backend.matrix)){
      backend.matrix[i,j] <- median(global.t[get.clus,j])
    }
  }
  hm_pal = heatmap_color_palette
  z <- backend.matrix
  col_seq <- seq(min(z),max(z), l = n <- 100)
  if(include_clusters[1]!='all') {
    backend.matrix <- backend.matrix[which(row.names(backend.matrix) %in% include_clusters),]
  }
  #color.map.fun = circlize::colorRamp2(seq(min(z),max(z), l = n <- 100), colorRampPalette(hm_pal)(n))
  color.map.fun = circlize::colorRamp2(col_seq, colorRampPalette(hm_pal)(n))
  # ncell <- rep(NA,times=length(unique(cluster_numbers)))
  ncell <- rep(NA,times=nrow(backend.matrix))
  # names(ncell) <- unique(cluster_numbers)[order(unique(cluster_numbers))]
  names(ncell) <- row.names(backend.matrix)
  for(i in 1:length(ncell)) {
    ncell[i] <- sum(cluster_numbers==names(ncell)[i])
  }
  pop.freq <- matrix(data=ncell,ncol=1)
  row.names(pop.freq) <- names(ncell)
  size_anno_nums <- round((pop.freq/sum(pop.freq))*100,2)
  ranno1 <- rowAnnotation(`Cluster\nSize`=anno_barplot(pop.freq,border=F,width=unit(1.75, "cm"),
                                                       axis_param=list(gp=gpar(fontsize=9)), axis = TRUE),
                          annotation_name_gp=gpar(fontsize=10,fontface="bold"), name = "Cluster\nSize")
  ranno2 <- rowAnnotation(frequency=anno_text(paste0(size_anno_nums,"%"),
                                              gp=gpar(fontsize=10,fontface="bold")))
  backend.matrix <- backend.matrix[order(row.names(backend.matrix)),]
  if(transpose_heatmap) {
    backend.matrix <- t(backend.matrix)
  }
  if(return_heatmap_data) {
    return(backend.matrix)
  }
  heatmap_output <- Heatmap(backend.matrix,col=color.map.fun,
                            row_names_side="left",
                            name="median\nscaled\nexpression",
                            heatmap_legend_param=list(at=c(0,0.2,0.4,0.6,0.8,1),legend_height=unit(3,"cm"),
                                                      grid_width=unit(0.6,"cm"),title_position="topleft",
                                                      labels_gp=gpar(fontsize=legend_text_size),
                                                      title_gp=gpar(fontsize=legend_text_size)),
                            row_names_gp=gpar(fontsize=row_text_size,fontface="bold"),
                            column_names_gp=gpar(fontsize=column_text_size,fontface="bold"),
                            rect_gp = gpar(lwd = heatmap_linewidth, col = "black"), border = "black",
                            cluster_columns = ifelse(cluster_col,TRUE,FALSE), cluster_rows = ifelse(cluster_row,TRUE,FALSE),
                            row_gap=unit(1,"mm"),column_gap=unit(1,"mm"),row_dend_gp=gpar(lwd=1.2),row_dend_width=unit(1,"cm"),
                            column_dend_gp = gpar(lwd=1.2), column_dend_height = unit(1,"cm")) +
    ranno1 + ranno2
  fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]] <- list(heatmap = heatmap_output,
                                                                heatmap_tile_data = backend.matrix,
                                                                population_size = pop.freq,
                                                                rep_used = ifelse(cordat,"with batch correction","without batch correction"))
  if(!'object_history' %in% names(fcs_join_obj)) {
    print("Consider running FCSimple::fcs_update() on the object.")
  } else {
    fcs_join_obj[['object_history']] <- append(fcs_join_obj[['object_history']],
                                               paste0(tolower(algorithm)," heatmap on ",ifelse(cordat,"corrected","uncorrected")," data: ",Sys.time()))
  }
  return(fcs_join_obj)
}

#' @title Save Cluster Heatmap to PDF
#'
#' @description
#'   Exports a ClusterHeatmap (generated by
#'   fcs_cluster_heatmap) to a PDF file. Filenames encode algorithm,
#'   correction status, and optional timestamp or suffix.
#'
#' @param fcs_join_obj
#'   A list with an element `<algorithm>_heatmap` as produced by
#'   FCSimple::fcs_cluster_heatmap().
#'
#' @param algorithm
#'   Character; name of the heatmap element to save (e.g. `"leiden"`).
#'
#' @param outdir
#'   Character; path to an existing directory. Defaults to `getwd()`.
#'
#' @param add_timestamp
#'   Logical; if `TRUE`, append a timestamp to the filename.
#'   Default `TRUE`.
#'
#' @param append_file_string
#'   Character or `NA`; if provided, a custom suffix to append to
#'   the PDF filename (before “.pdf”).
#'
#' @return
#'   Invisibly returns `NULL` after saving the PDF.
#'
#' @examples
#' \dontrun{
#'   # Assume `clustered` has a “leiden_heatmap” element
#'   FCSimple::fcs_plot_heatmap(
#'     clustered,
#'     algorithm = "leiden",
#'     outdir = "~/results",
#'     add_timestamp = FALSE,
#'     append_file_string = "v1"
#'   )
#' }
#'
#' @seealso
#'   FCSimple::fcs_cluster_heatmap
#'
#' @importFrom ggplot2 ggsave
#' @importFrom grid grid.grabExpr
#' @importFrom ComplexHeatmap draw
#' @export
fcs_plot_heatmap <- function(fcs_join_obj, algorithm, outdir = getwd(), add_timestamp = TRUE, append_file_string = NA)
{
  require(ggplot2)
  require(ComplexHeatmap)

  if(tolower(algorithm)=="dbscan") {
    fname <- paste0(outdir,"/",tolower(algorithm),"_cluster_heatmap_dbscan.pdf")
  } else {
    cordat <- ifelse(grepl(pattern = "with batch", x = fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][['rep_used']]),TRUE,FALSE)
    if(add_timestamp) {
      fname <- paste0(outdir,"/",tolower(algorithm),"_",ifelse(cordat,"cor","uncor"),"_cluster_heatmap_",strftime(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf")
    } else {
      fname <- paste0(outdir,"/",tolower(algorithm),"_",ifelse(cordat,"cor","uncor"),"_cluster_heatmap.pdf")
    }
  }
  if(!is.na(append_file_string)) {
    fname <- gsub(pattern = ".pdf$", replacement = paste0("_",append_file_string,".pdf"), x = fname)
  }
  ggsave(filename = fname,
         plot = grid::grid.grabExpr(draw(fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap"]])),
         device = "pdf", width = (ncol(fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]])*0.33)+3.25,
         height = (nrow(fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]])*0.3)+2.25,
         units = "in", dpi = 900)
}
