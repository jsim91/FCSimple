fcs_plot_distributions <- function(fcs_join_obj,
                                   separate_by = c("none", "date"))
{
  if(separate_by=="date"){
    if(!"run_date" %in% names(fcs_join_obj)){
      print("Unable to find run date. Using separate_by = 'none' instead.")
    } else {

    }
  }
}
