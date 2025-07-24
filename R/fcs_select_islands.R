#’ @title  
#’   Interactive DBSCAN-Based Island Selection in a 2D Embedding
#’
#’ @description  
#’   Performs density-based clustering (DBSCAN) on a UMAP or tSNE embedding,  
#’   renders a heatmap of cluster “islands” over selected lineage markers,  
#’   prompts the user to pick which islands to retain, and highlights those  
#’   cells on the original embedding. Finally returns a pared-down analysis  
#’   object containing only the selected events.
#’
#’ @param fcs_join_obj  
#’   A list returned by FCSimple::fcs_join() and FCSimple::fcs_reduce_dimensions(),  
#’   containing at minimum:  
#’   - `data`: numeric matrix of events × channels  
#’   - `source`: sample identifier vector  
#’   - `umap` or `tsne` element with `$coordinates` (cells × 2)  
#’
#’ @param reduction_cluster_annotate_algorithm  
#’   Character; clustering result to use when highlighting kept islands on the  
#’   embedding. One of `"leiden"`, `"flowsom"`, `"louvain"`, or `"phenograph"`.  
#’   Default uses the first value.
#’
#’ @param dbscan_reduction  
#’   Character; which reduction to run DBSCAN on. Either `"umap"` (default) or  
#’   `"tsne"`. Must match a named element in `fcs_join_obj`.
#’
#’ @param dbscan_minpts  
#’   Integer; minimum points per DBSCAN cluster (`minPts`). Default 50.
#’
#’ @param dbscan_eps  
#’   Numeric; epsilon neighborhood radius for DBSCAN. Default 0.5.
#’
#’ @param dbscan_lineage  
#’   Character vector of channel names to include in the island heatmap (e.g.  
#’   lineage markers). Default `c("170Yb_CD3","165Ho_CD19","145Nd_CD4","146Nd_CD8")`.  
#’   If `NULL`, uses all channels.
#’
#’ @param outdir  
#’   Character; path to an existing directory where intermediate heatmaps and  
#’   plots will be saved. Defaults to the working directory.
#’
#’ @details  
#’   1. Selects the specified embedding (`umap` or `tsne`) and if there are  
#’      more than 500,000 cells, randomly samples 500,000 for DBSCAN and uses  
#’      k-NN to classify the remainder.  
#’   2. Runs `dbscan::dbscan()` with `eps = dbscan_eps` and `minPts = dbscan_minpts`.  
#’   3. Builds a temporary FCSimple object containing only the DBSCAN clusters  
#’      and (optionally) the `dbscan_lineage` channels.  
#’   4. Calls `fcs_cluster_heatmap()` + `fcs_plot_heatmap()` on this object to  
#’      display median expression of each DBSCAN island.  
#’   5. Prompts the user via `readline()` to enter comma-separated cluster IDs  
#’      to keep.  
#’   6. Calls `fcs_plot_reduction()` with `internal_call = TRUE` to highlight  
#’      the selected events on the original embedding, colored by  
#’      `reduction_cluster_annotate_algorithm`.  
#’   7. Returns a pared object (list) with:  
#’      - `data`: the original expression matrix subset to kept rows  
#’      - `source`: subset of the `source` vector  
#’      - `run_date` (if present) subset to kept rows  
#’      - `subset_cells`: integer indices of kept events  
#’
#’ @return  
#’   A list with elements:  
#’   - `data`: numeric matrix of the selected events × channels  
#’   - `source`: character vector of sample IDs for selected events  
#’   - `run_date`: character vector of acquisition dates (if available)  
#’   - `subset_cells`: integer indices of the retained cells  
#’
#’ @examples  
#’ \dontrun{  
#’   # Perform UMAP reduction first  
#’   joined <- FCSimple::fcs_join(files)  
#’   reduced <- FCSimple::fcs_reduce_dimensions(joined, algorithm = "umap")  
#’  
#’   # Launch interactive island selector  
#’   sel <- FCSimple::fcs_select_islands(  
#’     reduced,  
#’     reduction_cluster_annotate_algorithm = "leiden",  
#’     dbscan_reduction   = "umap",  
#’     dbscan_minpts      = 40,  
#’     dbscan_eps         = 0.4,  
#’     dbscan_lineage     = c("CD3","CD19","CD4","CD8"),  
#’     outdir             = "~/my_analysis/plots"  
#’   )  
#’   # sel$data, sel$source, sel$run_date, sel$subset_cells  
#’ }  
#’  
#’ @seealso  
#’   dbscan::dbscan, FNN::knn,  
#’   FCSimple::fcs_cluster_heatmap, FCSimple::fcs_plot_heatmap,  
#’   FCSimple::fcs_plot_reduction  
#’  
#’ @importFrom dbscan dbscan  
#’ @importFrom FNN knn  
#’ @export
fcs_select_islands <- function(fcs_join_obj,
                               reduction_cluster_annotate_algorithm = c("leiden","flowsom","louvain","phenograph"),
                               dbscan_reduction = c("umap","tsne"),
                               dbscan_minpts = 50,
                               dbscan_eps = 0.5,
                               dbscan_lineage = c("170Yb_CD3","165Ho_CD19","145Nd_CD4","146Nd_CD8"),
                               outdir = getwd())
{
  require(dbscan)
  require(FNN)

  if(!tolower(dbscan_reduction) %in% names(fcs_join_obj)) {
    stop(paste0("Cannot find specified reduction ",dbscan_reduction," to cluster on. Have you run 'fcs_reduce_dimensions' yet?"))
  }
  if(dbscan_reduction == "umap") {
    dbscan_input <- fcs_join_obj[["umap"]][["coordinates"]]
  } else if(dbscan_reduction == "tsne") {
    dbscan_input <- fcs_join_obj[["tsne"]][["coordinates"]]
  } else {
    stop("error in argument 'dbscan_reduction': use either 'umap' or 'tsne'")
  }
  if(nrow(dbscan_input)>500000) {
    set.seed(123)
    sample_rows <- sample(1:nrow(dbscan_input),size=500000,replace=F)
    scan_out_sampled <- dbscan(x = dbscan_input[sample_rows,], eps = dbscan_eps, minPts = dbscan_minpts)
    time1 <- Sys.time()
    fnn_classify <- FNN::knn(train = dbscan_input[sample_rows,], test = dbscan_input[-sample_rows,],
                             cl = factor(scan_out_sampled$cluster), k = 5, prob = FALSE, algorithm = "kd_tree")
    # difftime(Sys.time(),time1,units="mins") # 6.5mins 0.5M -> 5.5M @ 5
    cluster_numbers <- rep(NA,times=nrow(dbscan_input))
    cluster_numbers[sample_rows] <- scan_out_sampled$cluster; cluster_numbers[-sample_rows] <- as.numeric(as.character(fnn_classify))
  } else {
    scan_out <- dbscan(x = dbscan_input, eps = dbscan_eps, minPts = dbscan_minpts)
    cluster_numbers <- scan_out$cluster
  }
  tmp_join_obj <- list(source = rep(NA,times=length(cluster_numbers)),
                       dbscan = list(clusters = cluster_numbers),
                       data = fcs_join_obj[["data"]])
  if(!is.null(dbscan_lineage)) {
    col_indices <- which(colnames(tmp_join_obj$data) %in% dbscan_lineage)
    if(length(col_indices)!=0) {
      tmp_join_obj[["data"]] <- tmp_join_obj[["data"]][,col_indices]
    }
  }
  calc_hm <- FCSimple::fcs_cluster_heatmap(fcs_join_obj = tmp_join_obj, algorithm = "dbscan")
  FCSimple::fcs_plot_heatmap(fcs_join_obj = calc_hm, algorithm = "dbscan", outdir = outdir)

  print("Using the dbscan heatmap here: ")
  print(paste0(paste0(gsub("/$","",outdir),"/dbscan_cluster_heatmap.pdf"),")"))
  user_input <- readline("which dbscan clusters should be kept? Enter integer values separated by commas:")
  keep_clus <- gsub(" ","",user_input)
  keep_clus <- as.numeric(strsplit(x = keep_clus, split = ",")[[1]])
  keep_clus <- unique(keep_clus)

  dbscan_keep_rows <- which(cluster_numbers %in% keep_clus)

  FCSimple::fcs_plot_reduction(fcs_join_obj = fcs_join_obj, algorithm = tolower(reduction_cluster_annotate_algorithm),
                               reduction = tolower(dbscan_reduction), point_alpha = 0.1, outdir = outdir,
                               internal_call = TRUE, keep_indices = dbscan_keep_rows)

  if(length(grep("DATE|date|Date",names(fcs_join_obj)))!=0) {
    fcs_obj_pared <- list(data = fcs_join_obj[["data"]][dbscan_keep_rows,],
                          source = fcs_join_obj[["source"]][dbscan_keep_rows],
                          run_date = fcs_join_obj[["run_date"]][dbscan_keep_rows],
                          subset_cells = dbscan_keep_rows)
  } else {
    fcs_obj_pared <- list(data = fcs_join_obj[["data"]][dbscan_keep_rows,],
                          source = fcs_join_obj[["source"]][dbscan_keep_rows],
                          subset_cells = dbscan_keep_rows)
  }
  return(fcs_obj_pared)
}
