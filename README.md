# FCSimple

Tools for the analysis of .fcs data sets. Pair with [FCView](https://github.com/jsim91/FCView) using `FCSimple::fcs_prepare_fcview_object` for an interactive analysis of FCSimple results.

## Installation

Install FCSimple using BiocManager to ensure all dependencies (including Bioconductor packages) are properly installed:

```r
# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install FCSimple with all dependencies
BiocManager::install("jsim91/FCSimple")
```

# Usage

This package requires the user to initiate a list object. If using csv file format:

```r
include_files <- list.files(path = "csv_file_directory", pattern = ".csv", full.names = TRUE)
my_object <- FCSimple::fcs_join(files = include_files)
```

or if using fcs format:

```r
include_files <- list.files(path = "fcs_file_directory", pattern = ".fcs", full.names = TRUE)
my_object <- FCSimple::fcs_join(files = include_files)
```

All subsequent functions take the object returned by `fcs_join()` as their first argument and return the updated object. The recommended minimum analysis order is shown below. Steps marked _optional_ may be skipped if not applicable to your data.

> **Note:** The pattern `my_object <- fcs_somefunction(my_object)` applies to all analysis and data-manipulation functions. **Visualisation functions are an exception** — they save figures to disk or return a plot object and should not be assigned back to `my_object`. Call them without assignment:
> ```r
> fcs_plot_distributions(fcs_join_obj = my_object)
> fcs_plot_reduction(fcs_join_obj = my_object)
> fcs_plot_reduction_density(fcs_join_obj = my_object)
> fcs_plot_reduction_difference(fcs_join_obj = my_object, compare_list = ..., color_list = ...)
> fcs_plot_overlays(fcs_join_obj = my_object)
> fcs_plot_trex(fcs_join_obj = my_object, ...)
> ```

```r
library(FCSimple)

# 1. Get the bundled example FCS files
include_files <- FCSimple::fcs_example_files()

# 2. Load data
my_object <- fcs_join(files = include_files)

# 3. Audit version info, validate metadata, and detect instrument type (optional)
my_object <- fcs_audit(fcs_join_obj = my_object)

# 4. Principal component analysis (optional) — used internally by batch correction
#    and as an optional input to clustering / dimension reduction
my_object <- fcs_pca(fcs_join_obj = my_object)

# 5. Batch correction — correct for acquisition-run effects (optional)
my_object <- fcs_batch_correction(fcs_join_obj = my_object)

# 6. Cluster cells
my_object <- fcs_cluster(fcs_join_obj = my_object, algorithm = "leiden")

# 7. Reduce dimensions for visualisation
my_object <- fcs_reduce_dimensions(fcs_join_obj = my_object, algorithm = "umap")

# 8. Calculate a cluster expression heatmap — required for the FCView heatmap tab
my_object <- fcs_cluster_heatmap(fcs_join_obj = my_object, algorithm = 'leiden')

# 9. Calculate cluster frequencies (proportion of each cluster per sample)
my_object <- fcs_calculate_abundance(fcs_join_obj = my_object,
                                     report_algorithm = 'leiden',
                                     report_as = "frequency")

# 10. Calculate cluster counts (number of cells per cluster per sample)
my_object <- fcs_calculate_abundance(fcs_join_obj = my_object,
                                     report_algorithm = 'leiden',
                                     report_as = "count")

# 11. Add sample-level metadata (clinical variables, group labels, etc.)
#     First column of added metadata must be 'patient_ID' and match
#     current metadata column 'patient_ID' to allow a proper merge.
my_metadata_df <- data.frame(patient_ID = my_object$metadata$patient_ID, 
                             visit = stringr::str_extract(my_object$metadata$patient_ID, 'v[1-3]'))
my_object <- fcs_add_metadata(fcs_join_obj = my_object, custom_metadata = my_metadata_df)

# 12. (Optional) Pre-annotate clusters with cell type labels before FCView upload
#     It is perfectly fine to only annotate a subset of the clusters at this point
#     The annotations made here are simply a proof of concept
my_object <- fcs_annotate_clusters(
  fcs_join_obj = my_object,
  clustering_algorithm = 'leiden',
  annotations = list(
    "T cell" = c(1, 3, 5),
    "B cell" = c(2, 4)
  )
)

# 13. Prepare and export for FCView
#     Note: no file extension on file_name; .RData added internally
#     Both saves to disk and returns to environment a formatted object
my_object_fcview <- fcs_prepare_fcview_object(
  fcs_join_obj = my_object,
  clustering_algorithm = 'leiden',
  output_dir = "path/to/output",
  file_name = "my_analysis"
)
```

See `?FCSimple::<function_name>` for full argument documentation on any step.

Once the object is exported, load it into **[FCView](https://github.com/jsim91/FCView)** for interactive visualisation, statistical testing, and annotation.

# Data Input Format

It's highly recommended to set transform parameters for all features in Flowjo then:

1) If feeding in .csv files to fcs_join(), you may prepare the .csv files by exporting the transformed data from Flowjo:

<img src="https://raw.githubusercontent.com/jsim91/FCSimple/main/vignette_outs/rmd_img/flowjo_csv_export_method.png" alt="FlowJo CSV export method">

2) If feeding in .fcs files to fcs_join(), you may prepare a workspace diagnostics file (as .txt file) for use with fcs_join() using the flowjo_diagnostics_file argument. To prepare the workspace diagnostics file, in Flowjo:

<img width="1588" alt="readme_flowjo_diagnostics" src="https://github.com/user-attachments/assets/88b137d4-9947-48e1-952f-035646626245">

3) If feeding in .fcs files without a diagnostics file, you can set `transform_per_channel = TRUE` in `fcs_join()` to launch an interactive Shiny app for inspecting and assigning transforms per channel. The app lets you preview density and 2D scatter plots under different transform settings (asinh, biexponential, hyperlog, or linear), apply settings to individual channels or all channels at once, and preview a sample UMAP before finalising. The chosen transform settings are stored on the returned object under `$app_transforms` for reproducibility:

```r
my_object <- fcs_join(files = include_files, transform_per_channel = TRUE)
my_object$app_transforms  # inspect the applied settings
```
<img width="1657" height="1487" alt="image" src="https://github.com/user-attachments/assets/168825e3-e195-45ea-b5ed-ce14cc4c6328" />

# Python Dependencies (optional)

The clustering and dimension reduction steps offer methods to run calculations through Python. No Python knowledge is required — functions call Python in the background and return results to R. The reticulate package is not required. Python scripts are located at `inst/python` and may be edited to alter default behaviour.

To take advantage of Python-supported methods, install the following:

- [Python](https://www.python.org/downloads/) — check "Add to PATH" during installation; pip is included
- [FlowKit](https://pypi.org/project/FlowKit/) — **required** for hyperlog transformation when using a FlowJo workspace diagnostics file with `fcs_join()` (not needed for the standard `transform = "hyperlog"` argument)
- [umap](https://github.com/lmcinnes/umap) — UMAP via Python
- [pynndescent](https://github.com/lmcinnes/pynndescent) — required for Python UMAP
- [leidenalg](https://github.com/vtraag/leidenalg) — Leiden clustering via Python
- [scipy](https://pypi.org/project/scipy/)
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [igraph](https://pypi.org/project/igraph/) — Louvain clustering via Python
- [openTSNE](https://github.com/pavlin-policar/openTSNE) — tSNE via Python

To install or verify all Python dependencies run:

```r
FCSimple::fcs_install_python_dependencies()
# or: ?FCSimple::fcs_install_python_dependencies
```
