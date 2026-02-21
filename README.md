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

```r
library(FCSimple)

# 1. Load data
my_object <- fcs_join(files = include_files)

# 2. Apply package / parameter updates to an existing object (optional)
my_object <- fcs_update(fcs_join_obj = my_object)

# 3. Principal component analysis — used internally by batch correction
#    and as an optional input to clustering / dimension reduction
my_object <- fcs_pca(fcs_join_obj = my_object)

# 4. Batch correction — correct for acquisition-run effects (optional)
my_object <- fcs_batch_correction(fcs_join_obj = my_object)

# 5. Cluster cells
my_object <- fcs_cluster(fcs_join_obj = my_object, algorithm = "leiden")

# 6. Reduce dimensions for visualisation
my_object <- fcs_reduce_dimensions(fcs_join_obj = my_object, algorithm = "umap")

# 7. Calculate a cluster expression heatmap (optional)
my_object <- fcs_cluster_heatmap(fcs_join_obj = my_object)

# 8. Calculate cluster frequencies (proportion of each cluster per sample)
my_object <- fcs_calculate_abundance(fcs_join_obj = my_object, report_as = "frequency")

# 9. Calculate cluster counts (number of cells per cluster per sample)
my_object <- fcs_calculate_abundance(fcs_join_obj = my_object, report_as = "count")

# 10. Attach sample-level metadata (clinical variables, group labels, etc.)
my_object <- fcs_add_metadata(fcs_join_obj = my_object, metadata = my_metadata_df)

# 11. (Optional) Pre-annotate clusters with cell type labels before FCView upload
my_object <- fcs_annotate_clusters(
  fcs_join_obj = my_object,
  annotations  = list(
    "T cell" = c(1, 3, 5),
    "B cell" = c(2, 4)
  )
)

# 12. Prepare and export for FCView
my_object_fcview <- fcs_prepare_fcview_object(
  fcs_join_obj = my_object,
  output_dir   = "path/to/output",
  file_name    = "my_analysis"
)
```

See `?FCSimple::<function_name>` for full argument documentation on any step.

# Data Input Format

It's highly recommended to set transform parameters for all features in Flowjo then:

1) If feeding in .csv files to fcs_join(), you may prepare the .csv files by exporting the transformed data from Flowjo:

<img src="https://raw.githubusercontent.com/jsim91/FCSimple/main/vignette_outs/rmd_img/flowjo_csv_export_method.png" alt="FlowJo CSV export method">

2) If feeding in .fcs files to fcs_join(), you may prepare a workspace diagnostics file (as .txt file) for use with fcs_join() using the flowjo_diagnostics_file argument. To prepare the workspace diagnostics file, in Flowjo:

<img width="1588" alt="readme_flowjo_diagnostics" src="https://github.com/user-attachments/assets/88b137d4-9947-48e1-952f-035646626245">

Tested in and developed with the Windows 10/11 operating system.

# Dependencies

## R packages

- [flowCore](https://www.bioconductor.org/packages/release/bioc/html/flowCore.html)
- [ncdfFlow](https://www.bioconductor.org/packages/release/bioc/html/ncdfFlow.html)
- [flowWorkspace](https://www.bioconductor.org/packages/release/bioc/html/flowWorkspace.html)
- [uwot](https://github.com/jlmelville/uwot)
- [Rtsne](https://github.com/jkrijthe/Rtsne)
- [FNN](https://cran.r-project.org/web/packages/FNN/index.html)
- [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
- [igraph](https://cran.r-project.org/web/packages/igraph/index.html)
- [FlowSOM](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)
- [Rphenograph](https://github.com/JinmiaoChenLab/Rphenograph)
- [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)
- [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
- [CATALYST](https://www.bioconductor.org/packages/release/bioc/html/CATALYST.html)
- [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [shadowtext](https://cran.r-project.org/web/packages/shadowtext/index.html)
- [future](https://cran.r-project.org/web/packages/future/index.html)
- [future.apply](https://cran.r-project.org/web/packages/future.apply/index.html)
- [RANN](https://cran.r-project.org/web/packages/RANN/RANN.pdf)
- [shiny](https://cran.r-project.org/web/packages/shiny/index.html)
- [shinythemes](https://cran.r-project.org/web/packages/shinythemes/index.html)
- [scales](https://cran.r-project.org/web/packages/scales/index.html)
- [cytoMEM](https://bioconductor.org/packages/release/bioc/html/cytoMEM.html)
- [dbscan](https://cran.r-project.org/web/packages/dbscan/index.html)
- [cyCombine](https://github.com/biosurf/cyCombine)

## Python (optional)

The clustering and dimension reduction steps offer methods to run calculations through Python. No Python knowledge is required — functions call Python in the background and return results to R. The reticulate package is not required. Python scripts are located at `inst/python` and may be edited to alter default behaviour.

To take advantage of Python-supported methods, install the following:

- [Python](https://www.python.org/downloads/) — check "Add to PATH" during installation; pip is included
- [FlowKit](https://pypi.org/project/FlowKit/) — **required** for hyperlog transformation via `fcs_join()`
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
