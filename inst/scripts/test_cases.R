fcs_test_file <- FCSimple::fcs_example_files()
fcs_obj <- FCSimple::fcs_join(files = fcs_test_file, instrument_type = 'cytof', use_descriptive_column_names = T,
                              transform_type = 'asinh', asinh_transform_cofactor = 5, transform_per_channel = F,
                              downsample_size = 100000)
fcs_obj$run_date <- 'batch1'
fcs_obj <- FCSimple::fcs_cluster(fcs_join_obj = fcs_obj, use_rep = 'data', language = 'python', algorithm = 'leiden',
                                 leiden_louvain_resolution = 1, adjacency_knn = 30, search_method = 'RANN',
                                 search_only = F, num_cores = 20)
fcs_obj <- FCSimple::fcs_reduce_dimensions(fcs_join_obj = fcs_obj, use_rep = 'data', algorithm = 'umap',
                                           language = 'python', umap_nn = 30, umap_min_dist = 0.1, nthread = 20)
fcs_obj <- FCSimple::fcs_calculate_abundance(fcs_join_obj = fcs_obj, report_algorithm = 'leiden', report_as = 'frequency')

umap_red <- FCSimple::fcs_plot_reduction(fcs_join_obj = fcs_obj, algorithm = 'leiden', reduction = 'umap',
                                         annotate_text_size = 5, return_plot = T, sample_equally = T)

fcs_obj <- FCSimple::fcs_cluster_heatmap(fcs_join_obj = fcs_obj, algorithm = 'leiden', override_correction = FALSE)

test_compare_list <- list('V1' = c('cytof_V1.fcs'),
                          'V2' = c('cytof_V2.fcs'),
                          'V3' = c('cytof_V3.fcs'))
test_color_list <- list('V1' = 'red',
                        'V2' = 'blue',
                        'V3' = 'purple')
test_comparisons_list <- list(c('V1','V2'),
                              c('V1', 'V3'),
                              c('V2','V3'))
fcs_obj <- FCSimple::fcs_test_clusters(fcs_join_obj = fcs_obj, compare_list = test_compare_list,
                                       color_list = test_color_list, comparisons = test_comparisons_list,
                                       denominator_cell_type = 'CD45+', algorithm = 'leiden', paired_test = F,
                                       p_text_size = 6)
fcs_obj$leiden$cluster_test_results[[3]]


fcs_join_obj <- fcs_obj
compare_list <- test_compare_list
color_list <- test_color_list
comparisons <- test_comparisons_list
denominator_cell_type <- 'CD45+'
x_order = NULL
abundance = NA
heatmap_matrix = NA
force_max = FALSE
algorithm = "leiden"
Rcolorbrewer_palette = "RdYlBu"
dot_size = 1
overlay_heatmap_numbers = TRUE
paired_test = FALSE
p_text_size = 5
paired_line_stroke = 0.1
paired_line_color = "black"
heatmap_fontsize = 8
relative_heights = c(0.76,0.24)
heatmap_parameters = 'all'
heatmap_clusters = "all"
test_method = 'sccomp'

input = in_list[[3]]
dplot_col = plot_cols
compare_these = my_compare
backmat = hm_tiles
use_palette = Rcolorbrewer_palette
size_of_dots = dot_size
cell_type_denom = denominator_cell_type
heatmap_overlay_values = overlay_heatmap_numbers
fm = force_max
abundance_alg = algorithm
pair_test = paired_test
xord = x_order
pts = p_text_size
pls = paired_line_stroke
plc = paired_line_color
relh = relative_heights
hmfs = heatmap_fontsize
tm = test_method


