% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_tests.R
\name{fcs_test_clusters}
\alias{fcs_test_clusters}
\title{Statistical Testing of Cluster Abundances}
\usage{
fcs_test_clusters(
  fcs_join_obj,
  compare_list,
  color_list,
  comparisons,
  denominator_cell_type,
  x_order = NULL,
  abundance = NA,
  heatmap_matrix = NA,
  force_max = FALSE,
  algorithm = c("leiden", "flowsom", "louvain", "phenograph", "git"),
  Rcolorbrewer_palette = "RdYlBu",
  dot_size = 1,
  overlay_heatmap_numbers = TRUE,
  paired_test = FALSE,
  p_text_size = 5,
  paired_line_stroke = 0.1,
  paired_line_color = "black",
  heatmap_fontsize = 8,
  relative_heights = c(0.76, 0.24),
  heatmap_parameters = "all",
  heatmap_clusters = "all"
)
}
\arguments{
\item{fcs_join_obj}{A list returned by FCSimple::fcs_join() and FCSimple::fcs_cluster(),
containing at minimum:
\itemize{
\item \code{fcs_join_obj[[ tolower(algorithm) ]][["abundance"]]}: sample × cluster matrix of abundances
\item (optionally) \code{fcs_join_obj[[ paste0(tolower(algorithm), "_heatmap") ]][["heatmap_tile_data"]]}
}}

\item{compare_list}{Named list of character vectors. Each element’s name defines a group
(e.g. “Control”, “Treatment”), and its value is the vector of sample
identifiers (matching row names of the abundance matrix) belonging to that group.}

\item{color_list}{Named character vector or list of colors for each group in \code{compare_list}.
Names must match those of \code{compare_list}.}

\item{comparisons}{List of length‐2 character vectors specifying which group pairs to test
(passed to ggpubr::stat_compare_means()). E.g. \code{list(c("Control","Case"))}.}

\item{denominator_cell_type}{Character; label used in the y‐axis title. The plot will read
\code{"\% of {denominator_cell_type}"}.}

\item{x_order}{Optional character vector giving the desired factor levels for the x‐axis.
If \code{NULL} (default), groups are taken in the order of \code{compare_list}.}

\item{abundance}{Optional numeric matrix of cluster abundances (samples × clusters). If
\code{NA} (default), pulled from
\code{fcs_join_obj[[ tolower(algorithm) ]][["abundance"]]}.}

\item{heatmap_matrix}{Optional numeric matrix of median expression values (clusters × markers).
If \code{NA} (default), pulled from the object’s existing heatmap tile data.}

\item{force_max}{Logical; if \code{TRUE}, and any frequency exceeds 95\%, the y‐axis max is set
to 100 (default \code{FALSE}).}

\item{algorithm}{Character; name of the clustering result to test (e.g. \code{"leiden"},
\code{"flowsom"}). Matches one of the elements in \code{fcs_join_obj}.}

\item{Rcolorbrewer_palette}{Character; name of an RColorBrewer palette of length ≥ 11 (default \code{"RdYlBu"}).}

\item{dot_size}{Numeric; size of dots in the dotplot (default \code{1}).}

\item{overlay_heatmap_numbers}{Logical; if \code{TRUE}, overlays numeric values on each heatmap cell (default \code{TRUE}).}

\item{paired_test}{Logical; if \code{TRUE}, adds paired Wilcoxon tests and connectors (default \code{FALSE}).}

\item{p_text_size}{Numeric; font size for p‐value labels in the distribution plot (default \code{5}).}

\item{paired_line_stroke}{Numeric; line width for paired connectors (default \code{0.1}).}

\item{paired_line_color}{Character; color for paired connectors (default \code{"black"}).}

\item{heatmap_fontsize}{Numeric; font size for numbers in the heatmap overlay (default \code{8}).}

\item{relative_heights}{Numeric vector of length 2 giving the vertical proportions of the
distribution plot vs. heatmap (default \code{c(0.76, 0.24)}).}

\item{heatmap_parameters}{Character vector of marker names to include in the heatmap (default \code{"all"}).}

\item{heatmap_clusters}{Character vector of cluster IDs to include in the heatmap (default \code{"all"}).}
}
\value{
The input \code{fcs_join_obj}, modified by adding:
\itemize{
\item \code{fcs_join_obj[[ tolower(algorithm) ]][["cluster_test_results"]]}:
a named list of ggarrange objects, one per cluster.
The function also prints a message indicating where the results are stored.
}
}
\description{
Compares cluster‐wise abundances between user‐defined sample groups and
visualizes results. For each cluster, the function generates a two‐panel
figure: a boxplot or dotplot with optional paired Wilcoxon tests, and a
companion heatmap of median marker expression. All cluster plots are
returned in a list within the input object for easy inspection or saving.
}
\details{
\enumerate{
\item Retrieve or use the provided \code{abundance} matrix and subset to the samples in \code{compare_list}.
\item For each cluster, construct a data.frame of frequencies with group labels.
\item Generate a distribution plot per cluster:
\itemize{
\item Boxplot or dotplot, colored by group
\item Optional paired or unpaired Wilcoxon tests via \code{ggpubr::stat_compare_means()}
\item Y-axis labeled \code{"\% of {denominator_cell_type}"}
\item If \code{force_max = TRUE}, force y-axis maximum at 100
}
\item Retrieve or compute a median-expression heatmap for the same clusters.
\item Assemble each cluster’s distribution plot + heatmap into a vertical two-panel figure using \code{ggpubr::ggarrange()}.
\item Store the resulting ggarrange list under
\code{fcs_join_obj[[ tolower(algorithm) ]][["cluster_test_results"]]}.
}
}
\examples{
\dontrun{
  files   <- list(ff1, ff2)
  joined  <- FCSimple::fcs_join(files)
  clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")

  tests <- FCSimple::fcs_test_clusters(
    clustered,
    compare_list = list(Healthy = c("S1","S2"), Diseased = c("S3","S4")),
    color_list   = c(Healthy = "steelblue", Diseased = "firebrick"),
    comparisons  = list(c("Healthy","Diseased")),
    denominator_cell_type = "T cells"
  )
  # View the figure for cluster "1":
  print(tests[["leiden"]][["cluster_test_results"]][["1"]])
}

}
\seealso{
FCSimple::fcs_cluster, FCSimple::fcs_calculate_abundance,
ggpubr::stat_compare_means, ComplexHeatmap::Heatmap
}
