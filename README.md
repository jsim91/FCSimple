todo:
- allow for downsampling to fit mclust then use in polygon test on full data to increase speed (cell gating)
- function that plots all gates set with hierachy (tree, parent, and current gate name included in title?) and writes to single multi-panel figure
- ensure tail gate works when tail < prominent peak
- test general gating functionality on multiple projects where fluorescence signature may differ
- resolve function conflicts
- ensure cells gate works on other flow projects where the distribution may be slightly different
note:
- set gate function should store a mask for both directions for later phenotype masking

# FCSimple

Tools for the analysis of .fcs data sets. Pair with [FCView](https://github.com/jsim91/FCView) using ```FCSimple::fcs_prepare_fcview_object``` for an interactive analysis of FCSimple results.

## Installation

Install FCSimple using BiocManager to ensure all dependencies (including Bioconductor packages) are properly installed:

```r
# Install BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install FCSimple with all dependencies
BiocManager::install("jsim91/FCSimple")
```

# Current dependencies include:

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

The dimension reduction and clustering steps do offer methods to run the calculations through Python, however no knowledge of the Python language is required. These functions will call Python in the background and results will be transferred to the R environment. The reticulate package is not required. Users may edit the included python scripts if they want to alter the default Python behavior. Scripts are located at: main/inst/python. To take advantage of Python-supported methods, these items are required:

- [Python](https://www.python.org/downloads/) # Check the "Add to Path" box during installation. The Python installation should come with pip.
- [FlowKit](https://pypi.org/project/FlowKit/) # *** necessary for hyperlog transformation of data when using fcs_join() ***
- [umap](https://github.com/lmcinnes/umap) # umap has its own set of dependencies, see link
- [pynndescent](https://github.com/lmcinnes/pynndescent) # Required to use this package's UMAP-by-Python functionality
- [leidenalg](https://github.com/vtraag/leidenalg) # leiden in Python
- [scipy](https://pypi.org/project/scipy/)
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [igraph](https://pypi.org/project/igraph/) # louvain in Python
- [openTSNE](https://github.com/pavlin-policar/openTSNE) # tSNE in Python


# Usage

This package requires the user to initiate a list object. If using csv file format:

```
include_files <- list.files(path = "csv_file_directory", pattern = ".csv", full.names = TRUE) # Point to the directory where the csv file(s) of interest are located. It's recommended to use full names.
my_object <- fcs_join(files = include_files) # This function takes the listed files as input and will initialize the list object that will store all subsequent analyses.
```

or if using fcs format:

```
include_files <- list.files(path = "fcs_file_directory", pattern = ".fcs", full.names = TRUE) # Point to the directory where the FCS file(s) of interest are located. It's recommended to use full names.
my_object <- fcs_join(files = include_files) # This function takes the listed files as input and will initialize the list object that will store all subsequent analyses.
```


Subsequent functions will take the object created with function

```
fcs_join()
```

as input. Analysis outputs will be joined to this object as new list entries. 

# Data Input Format

It's highly recommended to set transform parameters for all features in Flowjo then:

1) If feeding in .csv files to fcs_join(), you may prepare the .csv files by exporting the transformed data from Flowjo:

<img src="https://raw.githubusercontent.com/jsim91/FCSimple/main/vignette_outs/rmd_img/flowjo_csv_export_method.png" alt="FlowJo CSV export method">

2) If feeding in .fcs files to fcs_join(), you may prepare a workspace diagnostics file (as .txt file) for use with fcs_join() using the flowjo_diagnostics_file argument. To prepare the workspace diagnostics file, in Flowjo:

<img width="1588" alt="readme_flowjo_diagnostics" src="https://github.com/user-attachments/assets/88b137d4-9947-48e1-952f-035646626245">

Tested in and developed with the Windows 10/11 operating system.
