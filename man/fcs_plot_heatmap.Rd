% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_cluster_heatmap.R
\name{fcs_plot_heatmap}
\alias{fcs_plot_heatmap}
\title{Save Cluster Heatmap to PDF}
\usage{
fcs_plot_heatmap(
  fcs_join_obj,
  algorithm,
  outdir = getwd(),
  add_timestamp = TRUE,
  append_file_string = NA
)
}
\arguments{
\item{fcs_join_obj}{A list with an element \verb{<algorithm>_heatmap} as produced by
FCSimple::fcs_cluster_heatmap().}

\item{algorithm}{Character; name of the heatmap element to save (e.g. \code{"leiden"}).}

\item{outdir}{Character; path to an existing directory. Defaults to \code{getwd()}.}

\item{add_timestamp}{Logical; if \code{TRUE}, append a timestamp to the filename.
Default \code{TRUE}.}

\item{append_file_string}{Character or \code{NA}; if provided, a custom suffix to append to
the PDF filename (before “.pdf”).}
}
\value{
Invisibly returns \code{NULL} after saving the PDF.
}
\description{
Exports a ClusterHeatmap (generated by
fcs_cluster_heatmap) to a PDF file. Filenames encode algorithm,
correction status, and optional timestamp or suffix.
}
\examples{
\dontrun{
  # Assume `clustered` has a “leiden_heatmap” element
  FCSimple::fcs_plot_heatmap(
    clustered,
    algorithm = "leiden",
    outdir = "~/results",
    add_timestamp = FALSE,
    append_file_string = "v1"
  )
}

}
\seealso{
FCSimple::fcs_cluster_heatmap
}
