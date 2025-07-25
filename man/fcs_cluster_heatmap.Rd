% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_cluster_heatmap.R
\name{fcs_cluster_heatmap}
\alias{fcs_cluster_heatmap}
\title{Generate Cluster Heatmap for Flow Cytometry Data}
\usage{
fcs_cluster_heatmap(
  fcs_join_obj,
  algorithm,
  include_parameters = "all",
  include_clusters = "all",
  heatmap_color_palette = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
  transpose_heatmap = FALSE,
  cluster_row = TRUE,
  cluster_col = TRUE,
  override_correction = TRUE,
  return_heatmap_data = FALSE,
  heatmap_linewidth = 0.5
)
}
\arguments{
\item{fcs_join_obj}{A list returned by FCSimple::fcs_join() and FCSimple::fcs_cluster(),
containing at least:
\itemize{
\item \code{data}: numeric matrix of events × channels
\item \verb{<algorithm>$clusters}: vector of cluster IDs per event
\item optionally \code{batch_correction$data} if correction was applied
}}

\item{algorithm}{Character; name of clustering results in \code{fcs_join_obj} to use
(e.g. \code{"leiden"}, \code{"flowsom"}, etc.).}

\item{include_parameters}{Character vector of channel names to include (default \code{"all"}).
If \code{"all"}, uses all columns of the expression matrix.}

\item{include_clusters}{Character vector of cluster IDs to include (default \code{"all"}).
If \code{"all"}, includes every cluster.}

\item{heatmap_color_palette}{Character vector of colors (length ≥ 2) for the heatmap palette.
Default uses a reversed “RdYlBu” from RColorBrewer.}

\item{transpose_heatmap}{Logical; if \code{TRUE}, transpose the heatmap matrix before plotting.
Default \code{FALSE}.}

\item{cluster_row}{Logical; if \code{TRUE}, apply hierarchical clustering to rows.
Default \code{TRUE}.}

\item{cluster_col}{Logical; if \code{TRUE}, apply hierarchical clustering to columns.
Default \code{TRUE}.}

\item{override_correction}{Logical; if \code{TRUE}, always use raw data even if batch correction
is present. Default \code{TRUE}.}

\item{return_heatmap_data}{Logical; if \code{TRUE}, return the matrix of median‐scaled values
invisibly instead of adding the heatmap to \code{fcs_join_obj}.
Default \code{FALSE}.}

\item{heatmap_linewidth}{Numeric; border line width for heatmap cells. Default \code{0.5}.}
}
\value{
Invisibly returns the updated \code{fcs_join_obj} with a new element
\verb{<algorithm>_heatmap} as described above. If
\code{return_heatmap_data = TRUE}, returns only the heatmap matrix.
}
\description{
Computes and visualizes median‐scaled expression of each cluster across
parameters using ComplexHeatmap. Selects raw or batch‐corrected data,
aggregates by cluster, and stores both the heatmap object and tile data.
}
\details{
\itemize{
\item Chooses raw or batch‐corrected data based on \code{override_correction} and
presence of \code{fcs_join_obj$batch_correction$data}.
\item Scales expression values to \link{0,1} per channel using
CATALYST:::.scale_exprs.
\item Computes median expression for each cluster × parameter.
\item Builds a ComplexHeatmap object with cluster‐size annotations.
\item Stores the result under
\code{fcs_join_obj[[paste0(tolower(algorithm), "_heatmap")]]}:
\itemize{
\item \code{heatmap}: the Heatmap object
\item \code{heatmap_tile_data}: the numeric matrix used
\item \code{population_size}: cluster event counts
\item \code{rep_used}: “with batch correction” or “without batch correction”
}
\item Appends a timestamped entry to \code{object_history}.
}
}
\examples{
\dontrun{
  joined <- FCSimple::fcs_join(list(ff1, ff2))
  clustered <- FCSimple::fcs_cluster(joined, algorithm = "leiden")

  # Generate and store heatmap object
  out <- FCSimple::fcs_cluster_heatmap(
    clustered,
    algorithm = "leiden",
    include_parameters = c("CD3","CD4","CD8"),
    override_correction = FALSE
  )

  # Just get the tile matrix
  mat <- FCSimple::fcs_cluster_heatmap(
    clustered,
    algorithm = "leiden",
    return_heatmap_data = TRUE
  )
}

}
\seealso{
FCSimple::fcs_cluster, FCSimple::fcs_plot_heatmap
}
