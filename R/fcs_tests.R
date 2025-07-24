#’ @title  
#’   Statistical Testing of Cluster Abundances  
#’
#’ @description  
#’   Compares cluster‐wise abundances between user‐defined sample groups and  
#’   visualizes results. For each cluster, the function generates a two‐panel  
#’   figure: a boxplot or dotplot with optional paired Wilcoxon tests, and a  
#’   companion heatmap of median marker expression. All cluster plots are  
#’   returned in a list within the input object for easy inspection or saving.
#’
#’ @param fcs_join_obj  
#’   A list returned by FCSimple::fcs_join() and FCSimple::fcs_cluster(),  
#’   containing at minimum:  
#’   - `fcs_join_obj[[ tolower(algorithm) ]][["abundance"]]`: sample × cluster matrix of abundances  
#’   - (optionally) `fcs_join_obj[[ paste0(tolower(algorithm), "_heatmap") ]][["heatmap_tile_data"]]`  
#’
#’ @param compare_list  
#’   Named list of character vectors. Each element’s name defines a group  
#’   (e.g. “Control”, “Treatment”), and its value is the vector of sample  
#’   identifiers (matching row names of the abundance matrix) belonging to that group.
#’
#’ @param color_list  
#’   Named character vector or list of colors for each group in `compare_list`.  
#’   Names must match those of `compare_list`.
#’
#’ @param comparisons  
#’   List of length‐2 character vectors specifying which group pairs to test  
#’   (passed to ggpubr::stat_compare_means()). E.g. `list(c("Control","Case"))`.
#’
#’ @param denominator_cell_type  
#’   Character; label used in the y‐axis title. The plot will read  
#’   `"% of {denominator_cell_type}"`.
#’
#’ @param x_order  
#’   Optional character vector giving the desired factor levels for the x‐axis.  
#’   If `NULL` (default), groups are taken in the order of `compare_list`.
#’
#’ @param abundance  
#’   Optional numeric matrix of cluster abundances (samples × clusters). If  
#’   `NA` (default), pulled from  
#’   `fcs_join_obj[[ tolower(algorithm) ]][["abundance"]]`.
#’
#’ @param heatmap_matrix  
#’   Optional numeric matrix of median expression values (clusters × markers).  
#’   If `NA` (default), pulled from the object’s existing heatmap tile data.
#’
#’ @param force_max  
#’   Logical; if `TRUE`, and any frequency exceeds 95%, the y‐axis max is set  
#’   to 100 (default `FALSE`).
#’
#’ @param algorithm  
#’   Character; name of the clustering result to test (e.g. `"leiden"`,  
#’   `"flowsom"`). Matches one of the elements in `fcs_join_obj`.
#’
#’ @param Rcolorbrewer_palette  
#’   Character; name of an RColorBrewer palette of length ≥ 11 (default `"RdYlBu"`).
#’
#’ @param dot_size  
#’   Numeric; size of dots in the dotplot (default `1`).
#’
#’ @param overlay_heatmap_numbers  
#’   Logical; if `TRUE`, overlays numeric values on each heatmap cell (default `TRUE`).
#’
#’ @param paired_test  
#’   Logical; if `TRUE`, adds paired Wilcoxon tests and connectors (default `FALSE`).
#’
#’ @param p_text_size  
#’   Numeric; font size for p‐value labels in the distribution plot (default `5`).
#’
#’ @param paired_line_stroke  
#’   Numeric; line width for paired connectors (default `0.1`).
#’
#’ @param paired_line_color  
#’   Character; color for paired connectors (default `"black"`).
#’
#’ @param heatmap_fontsize  
#’   Numeric; font size for numbers in the heatmap overlay (default `8`).
#’
#’ @param relative_heights  
#’   Numeric vector of length 2 giving the vertical proportions of the  
#’   distribution plot vs. heatmap (default `c(0.76, 0.24)`).
#’
#’ @param heatmap_parameters  
#’   Character vector of marker names to include in the heatmap (default `"all"`).
#’
#’ @param heatmap_clusters  
#’   Character vector of cluster IDs to include in the heatmap (default `"all"`).
#’
#’ @details  
#’   1. Retrieves or uses the provided `abundance` matrix and subsets it to  
#’      the samples in `compare_list`.  
#’   2. Constructs, for each cluster, a data frame of frequencies with group labels.  
#’   3. Generates a distribution plot per cluster:  
#’      - Boxplot or dotplot, colored by group.  
#’      - Optional paired or unpaired Wilcoxon tests via `stat_compare_means()`.  
#’      - Y‐axis labeled with the specified cell‐type denominator.  
#’      - Optional y‐axis max at 100 when `force_max = TRUE`.  
#’   4. Retrieves or computes a median‐expression heatmap for the same clusters.  
#’   5. Assembles each cluster’s distribution plot and heatmap into a vertical  
#’      two‐panel figure using `ggpubr::ggarrange()`.  
#’   6. Stores the resulting list of ggarrange objects in  
#’      `fcs_join_obj[[ tolower(algorithm) ]][["cluster_test_results"]]`.
#’
#’ @return  
#’   The input `fcs_join_obj`, modified by adding:  
#’   - `fcs_join_obj[[ tolower(algorithm) ]][["cluster_test_results"]]`:  
#’     a named list of ggarrange objects, one per cluster.  
#’   The function also prints a message indicating where the results are stored.
#’
#’ @examples  
#’ \dontrun{  
#’   files   <- list(ff1, ff2)  
#’   joined  <- FCSimple::fcs_join(files)  
#’   clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")  
#’   
#’   tests <- FCSimple::fcs_test_clusters(  
#’     clustered,  
#’     compare_list = list(Healthy = c("S1","S2"), Diseased = c("S3","S4")),  
#’     color_list   = c(Healthy = "steelblue", Diseased = "firebrick"),  
#’     comparisons  = list(c("Healthy","Diseased")),  
#’     denominator_cell_type = "T cells"  
#’   )  
#’   # View the figure for cluster "1":  
#’   print(tests[["leiden"]][["cluster_test_results"]][["1"]])  
#’ }  
#’
#’ @seealso  
#’   FCSimple::fcs_cluster, FCSimple::fcs_calculate_abundance,  
#’   ggpubr::stat_compare_means, ComplexHeatmap::Heatmap  
#’
#’ @importFrom ggplot2 ggplot aes geom_boxplot geom_dotplot geom_line  
#’   scale_fill_manual scale_color_manual scale_y_continuous ylim labs  
#’   theme_minimal theme element_text  
#’ @importFrom ggpubr stat_compare_means ggarrange  
#’ @importFrom ComplexHeatmap Heatmap  
#’ @importFrom circlize colorRamp2  
#’ @importFrom grid text gpar grabExpr
#’ @export
fcs_test_clusters <- function(fcs_join_obj, compare_list, color_list, comparisons, denominator_cell_type,
                              x_order = NULL, abundance = NA, heatmap_matrix = NA, force_max = FALSE,
                              algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                              Rcolorbrewer_palette = "RdYlBu", # must be a colorbrewer palette that's 11 long such as Spectral or RdYlBu
                              dot_size = 1, overlay_heatmap_numbers = TRUE, paired_test = FALSE, 
                              p_text_size = 5, paired_line_stroke = 0.1, paired_line_color = "black", 
                              heatmap_fontsize = 8, relative_heights = c(0.76,0.24), heatmap_parameters = 'all', 
                              heatmap_clusters = "all")
{
  require(ggplot2)
  require(ggpubr)
  # require(ggtext)
  require(ggplot2)
  require(ggpubr)
  require(formattable)
  require(grid)
  require(ComplexHeatmap)
  require(circlize)

  # testing

  # end testing argument assignments

  if(is.na(abundance[1])) {
    abundance <- fcs_join_obj[[tolower(algorithm)]][["abundance"]]
  }
  rm_row <- which(!row.names(abundance) %in% unlist(compare_list))
  if(length(rm_row)!=0) {
    abundance <- abundance[-rm_row,]
  }
  for(i in 1:length(compare_list)) {
    try_rm <- which(!compare_list[[i]] %in% row.names(abundance))
    if(length(try_rm)!=0) {
      compare_list[[i]] <- compare_list[[i]][-try_rm]
    }
  }

  rbind_list <- vector("list", length = length(compare_list))
  for(i in 1:length(compare_list)) {
    rbind_list[[i]] <- t(t(abundance)[,compare_list[[i]]])
    rbind_list[[i]] <- cbind(rbind_list[[i]],data.frame(group = as.character(1:nrow(rbind_list[[i]]))))
  }

  group_assignment <- c()
  for(i in 1:length(compare_list)) {
    group_assignment <- append(group_assignment,rep(names(compare_list)[i],length(compare_list[[i]])))
  }

  abundance <- do.call(rbind,rbind_list)
  group_vals <- abundance[,"group"]
  abundance <- abundance[,-which(colnames(abundance)=="group")]

  cluster_list <- vector("list", length = ncol(abundance))
  names(cluster_list) <- colnames(abundance)
  for(i in 1:length(cluster_list)) {
    tmp_df <- data.frame(cluster = rep(colnames(abundance)[i],nrow(abundance)),
                         pid = row.names(abundance), group = group_vals,
                         frequency = abundance[,i])
    tmp_df$compare_group <- group_assignment
    cluster_list[[i]] <- tmp_df
  }

  # test_out <- vector("list", length = length(cluster_list)); names(test_out) <- names(cluster_list)
  plot_cols <- unlist(color_list)
  my_compare <- comparisons

  if(is.na(heatmap_matrix[1])) {
    if(heatmap_parameters[1]=='all') {
      hm_tiles <- fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]]
    } else {
      fcs_join_obj <- FCSimple::fcs_cluster_heatmap(fcs_join_obj = fcs_join_obj, algorithm = algorithm, include_parameters = heatmap_parameters, 
                                                    heatmap_color_palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu")), transpose_heatmap = FALSE, 
                                                    cluster_row = TRUE, cluster_col = TRUE, override_correction = TRUE, return_heatmap_data = FALSE, 
                                                    heatmap_linewidth = 0.5, include_clusters = heatmap_clusters)
      hm_tiles <- fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]]
    }
  } else {
    hm_tiles <- heatmap_matrix
  }
  in_list <- vector("list", length = length(cluster_list)); names(in_list) <- names(cluster_list)
  for(i in 1:length(cluster_list)) {
    in_list[[i]] <- list(cluster_list[[i]],t(as.matrix(hm_tiles[which(row.names(hm_tiles)==names(in_list)[i]),])))
  }

  test_plot <- function(input, dplot_col = plot_cols, compare_these = my_compare,
                        backmat = hm_tiles, use_palette = Rcolorbrewer_palette, 
                        size_of_dots = dot_size, cell_type_denom = denominator_cell_type,
                        heatmap_overlay_values = overlay_heatmap_numbers, fm = force_max,
                        abundance_alg = algorithm, pair_test = paired_test, xord = x_order, 
                        pts = p_text_size, pls = paired_line_stroke, plc = paired_line_color, 
                        relh = relative_heights, hmfs = heatmap_fontsize) {

    # testing #
    # abundance_alg <- algorithm
    # heatmap_overlay_values = overlay_heatmap_numbers
    # size_of_dots = dot_size
    # cell_type_denom = denominator_cell_type
    # input <- in_list[[1]]
    # compare_these = my_compare; dplot_col = plot_cols
    # use_palette = Rcolorbrewer_palette
    # pair_test = paired_test

    plot_input <- input[[1]]
    hm_input <- input[[2]]

    if(!is.null(xord)) {
      plot_input$compare_group <- factor(plot_input$compare_group, levels = xord)
    } else {
      plot_input$compare_group <- factor(plot_input$compare_group)
    }

    hm_pal <- rev(RColorBrewer::brewer.pal(11, use_palette))
    color.map.fun = circlize::colorRamp2(seq(0,1, l = n <- 100), colorRampPalette(hm_pal)(n))

    cluster_col_ind <- grep(pattern = "^cluster", x = colnames(plot_input))
    capture_cluster <- plot_input[1,cluster_col_ind]

    plt <- ggplot(data = plot_input, mapping = aes(x = compare_group, y = frequency, color = compare_group)) +
      geom_boxplot(fill = "#bfbfbf", lwd = 0.5, alpha = 0.4, width = 0.4)
    if(pair_test) {
      plt <- plt + geom_line(mapping = aes(group = group), color = "black", linewidth = pls) +
        geom_point(mapping = aes(fill = compare_group), pch = 21, size = size_of_dots*3, color = plc) +
        stat_compare_means(method = "wilcox", comparisons = compare_these, size = pts, paired = TRUE)
    } else {
      plt <- plt + geom_dotplot(data = plot_input, aes(fill = compare_group), color = plc, stroke = 0.5,
                                binaxis = "y", stackdir = "center", position = "dodge", binpositions="all",
                                dotsize = size_of_dots) +
        stat_compare_means(method = "wilcox", comparisons = compare_these, size = pts, paired = FALSE)
    }
    plt <- plt + scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(y = paste0("% of ",cell_type_denom), title = paste0("cluster ",capture_cluster)) +
      scale_fill_manual(values = dplot_col) +
      scale_color_manual(values = dplot_col) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14))
    if(fm) {
      if(max(plot_input$frequency)>95) {
        plt <- plt + ylim(floor(min(plot_input$frequency)), 100)
      }
    }

    colnames(hm_input) <- gsub("^.+_","",colnames(hm_input))
    if(heatmap_overlay_values) {
      hm_sub <- ComplexHeatmap::Heatmap(matrix = hm_input, col = color.map.fun, cluster_columns = FALSE,
                                        show_heatmap_legend = FALSE,
                                        cell_fun = function(j, i, x, y, width, height, fill) {
                                          grid.text(sprintf("%.2f", hm_input[i, j]), x, y, gp = gpar(fontsize = hmfs))
                                        })
    } else {
      hm_sub <- ComplexHeatmap::Heatmap(matrix = hm_input, col = color.map.fun, cluster_columns = FALSE,
                                        show_heatmap_legend = FALSE)
    }

    arr_plot <- ggpubr::ggarrange(plotlist = list(plt, grid.grabExpr(draw(hm_sub))), nrow = 2, ncol = 1,
                                  heights = relh)
    return(arr_plot)
  }

  out_plots <- lapply(X = in_list, FUN = test_plot)

  fcs_join_obj[[algorithm]][["cluster_test_results"]] <- out_plots

  print(paste0("test results saved in object under: fcs_join_obj[['",algorithm,"']][['cluster_test_results']]"))
  return(fcs_join_obj)
}
