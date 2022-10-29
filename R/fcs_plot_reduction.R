fcs_plot_reduction <- function(fcs_join_obj, algorithm, reduction)
{
  reduction_coords <- fcs_join_obj[[tolower(reduction)]][["coordinates"]]
  cluster_numbers <- fcs_join_obj[[tolower(algorithm)]][["clusters"]]
}
