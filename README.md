# FCSimple

Tools for the analysis of multidimensional FCS-based data sets. This is in active development. Bugs are likely and documentation will come later. Tasks that are currently supported include: 

```
fcs_join # Use for concatenation of multiple FCS files (can also take a single FCS file) and formatting of package-wide list object
fcs_reduce_dimensions # Use after running fcs_join() to produce low dimension representations of your data set with the UMAP or tSNE algorithms.
fcs_cluster # Use after running fcs_join() to partition the data into discrete clusters using the Leiden, Louvain, FlowSOM, or Phenograph algorithms.
fcs_cluster_abundance # Use after running fcs_cluster() to calculate cluster abundance, either as fraction (0-1) or frequency (0-100)
fcs_report_abundance # Use after running fcs_calculate_abundance() to report calculated cluster abundance values

# more to come
```


# Package Install

Enter the following commands within an R session

```
library(devtools)
install_github("jsim91/FCSimple")
```


# Current dependencies include:

- [flowCore](https://www.bioconductor.org/packages/release/bioc/html/flowCore.html)
- [ncdfFLow](https://www.bioconductor.org/packages/release/bioc/html/ncdfFlow.html)
- [flowWorkspace](https://www.bioconductor.org/packages/release/bioc/html/flowWorkspace.html)
- [uwot](https://github.com/jlmelville/uwot)
- [Rtsne](https://github.com/jkrijthe/Rtsne)
- [FNN](https://cran.r-project.org/web/packages/FNN/index.html)
- [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
- [leiden](https://github.com/TomKellyGenetics/leiden)
- [igraph](https://cran.r-project.org/web/packages/igraph/index.html)
- [FlowSOM](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)
- [Rphenograph](https://github.com/JinmiaoChenLab/Rphenograph)


The dimension reduction and clustering steps do offer methods to run the calculations through Python, however no knowledge of the Python language is required. These functions will call Python in the background and results will be transferred to the R environment. The reticulate package is not required. Users may edit the included python scripts if they want to alter the default Python behavior. To take davantage of Python-supported methods, these Python dependencies are required:

- [Python](https://www.python.org/downloads/) # Check the "Add to Path" box during installation. The Python installation should come with pip.
- [umap](https://github.com/lmcinnes/umap) # umap has its own set of dependencies, see link
- [pynndescent](https://github.com/lmcinnes/pynndescent) # Required to use this package's UMAP-by-Python functionality
- [leidenalg](https://github.com/vtraag/leidenalg)
- [scipy](https://pypi.org/project/scipy/)
- [numpy](https://pypi.org/project/numpy/)
- [pandas](https://pypi.org/project/pandas/)
- [igraph](https://pypi.org/project/igraph/) # Louvain is run using igraph


# Usage

This package requires the user to initiate a list object like this:

```
include_files <- list.files(path = "fcs_file_directory", full.names = TRUE) # Point to the directory where the FCS file(s) of interest are located. It's recommended to use full names.
my_object <- fcs_join(files = include_files) # This function takes the listed files as input and will initialize the list object that will store all subsequent analyses.
```

Subsequent functions will take the object created with function

```
fcs_join()
```

as input. Analysis outputs will be joined to this object as new list entries.
