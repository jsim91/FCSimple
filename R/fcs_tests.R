fcs_test_clusters <- function(fcs_join_obj, compare_list, color_list, comparisons,
                              algorithm = c("leiden","flowsom","louvain","phenograph","git"),
                              Rcolorbrewer_palette = "RdYlBu", # must be a colorbrewer palette that's 11 long such as Spectral or RdYlBu
                              dot_size = 1, denominator_cell_type, overlay_heatmap_numbers = TRUE)
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

  # testing argument assignments
  # overlay_heatmap_numbers = TRUE
  # denominator_cell_type = "lymphocytes"
  # Rcolorbrewer_palette = "RdYlBu"
  # dot_size = 1
  # fcs_join_obj <- readRDS("J:/Mashayekhi/flow/AHA_liraglutide_scenith/lymphocytes/AHA_liraglutide_lymphocytes_fcsimple_obj_clustered.rds")
  # algorithm <- "leiden"
  # abundance <- fcs_join_obj[[algorithm]][["abundance"]]
  # compare_list <- list('2DG' = row.names(abundance)[grep(pattern = "2DG", x = row.names(abundance))],
  #                      'DGO' = row.names(abundance)[grep(pattern = "DGO", x = row.names(abundance))],
  #                      'No Inhib' = row.names(abundance)[grep(pattern = "No Inhib", x = row.names(abundance))])
  # color_list <- list('2DG' = "red", 'DGO' = "blue", 'No Inhib' = "green")
  # comparisons <- list(c('2DG','DGO'),c('DGO','No Inhib'))
  # end testing argument assignments

  rm_row <- which(!row.names(abundance) %in% unlist(compare_list))
  if(length(rm_row)!=0) {
    abundance <- abundance[-rm_row,]
  }

  group_assignment <- rep(NA,nrow(abundance))
  for(i in 1:length(compare_list)) {
    group_assignment[which(row.names(abundance) %in% compare_list[[i]])] <- names(compare_list)[i]
  }

  cluster_list <- vector("list", length = ncol(abundance))
  names(cluster_list) <- colnames(abundance)
  for(i in 1:length(cluster_list)) {
    tmp_df <- data.frame(PID = row.names(abundance), clus = abundance[,names(cluster_list)[i]])
    colnames(tmp_df)[2] <- paste0("cluster_",names(cluster_list)[i])
    tmp_df$compare_group <- group_assignment
    cluster_list[[i]] <- tmp_df
  }

  test_out <- vector("list", length = length(cluster_list)); names(test_out) <- names(cluster_list)
  plot_cols <- unlist(color_list)
  my_compare <- comparisons

  hm_tiles <- fcs_join_obj[[paste0(tolower(algorithm),"_heatmap")]][["heatmap_tile_data"]]
  in_list <- vector("list", length = length(cluster_list)); names(in_list) <- names(cluster_list)
  for(i in 1:length(cluster_list)) {
    in_list[[i]] <- list(cluster_list[[i]],t(as.matrix(hm_tiles[which(row.names(hm_tiles)==names(in_list)[i]),])))
  }

  test_plot <- function(input, dplot_col = plot_cols, compare_these = my_compare,
                        backmat = hm_tiles, use_palette = Rcolorbrewer_palette,
                        size_of_dots = dot_size, cell_type_denom = denominator_cell_type,
                        heatmap_overlay_values = overlay_heatmap_numbers,
                        abundance_alg = algorithm) {

    abundance_alg <- algorithm
    heatmap_overlay_values = overlay_heatmap_numbers
    size_of_dots = dot_size
    cell_type_denom = denominator_cell_type
    input <- in_list[[1]]
    compare_these = my_compare; dplot_col = plot_cols
    use_palette = Rcolorbrewer_palette

    plot_input <- input[[1]]
    hm_input <- input[[2]]

    hm_pal <- rev(RColorBrewer::brewer.pal(11, use_palette))
    color.map.fun = circlize::colorRamp2(seq(0,1, l = n <- 100), colorRampPalette(hm_pal)(n))

    cluster_col_ind <- grep(pattern = "^cluster", x = colnames(plot_input))
    capture_cluster <- colnames(plot_input)[cluster_col_ind]
    colnames(plot_input)[cluster_col_ind] <- "placeholder"

    plt <- ggplot(data = plot_input, mapping = aes(x = group_assignment, y = placeholder, color = group_assignment)) +
      geom_boxplot(fill = "#bfbfbf", lwd = 0.5, alpha = 0.4, width = 0.4) +
      geom_dotplot(data = plot_input, aes(fill = group_assignment), color = "black", stroke = 1,
                   binaxis = "y", stackdir = "center", position = "dodge", binpositions="all",
                   dotsize = size_of_dots) +
      stat_compare_means(method = "wilcox", comparisons = compare_these, size = 4.5, paired = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      labs(y = paste0("% of ",cell_type_denom), title = gsub("_"," ",capture_cluster)) +
      scale_fill_manual(values = dplot_col) +
      scale_color_manual(values = dplot_col) +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            legend.position = "none",
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.text.x = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14))

    colnames(hm_input) <- gsub("^.+_","",colnames(hm_input))
    if(heatmap_overlay_values) {
      hm_sub <- ComplexHeatmap::Heatmap(matrix = hm_input, col = color.map.fun, cluster_columns = FALSE,
                                        show_heatmap_legend = FALSE,
                                        cell_fun = function(j, i, x, y, width, height, fill) {
                                          grid.text(sprintf("%.2f", hm_input[i, j]), x, y, gp = gpar(fontsize = 8))
                                        })
    } else {
      hm_sub <- ComplexHeatmap::Heatmap(matrix = hm_input, col = color.map.fun, cluster_columns = FALSE,
                                        show_heatmap_legend = FALSE)
    }

    arr_plot <- ggpubr::ggarrange(plotlist = list(plt, grid.grabExpr(draw(hm_sub))), nrow = 2, ncol = 1,
                                  heights = c(0.76,0.24))
    return(arr_plot)
  }

  out_plots <- lapply(X = in_list, FUN = test_plot)

  pdf(file = paste0(abundance_alg,"_cluster_statistics.pdf"), width = 10, height = 5)
  lapply(X = out_plots, FUN = function(x) x)
  dev.off()
}
