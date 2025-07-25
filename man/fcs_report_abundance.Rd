% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_cluster_abundance.R
\name{fcs_report_abundance}
\alias{fcs_report_abundance}
\title{Report Cluster Abundance to CSV File}
\usage{
fcs_report_abundance(
  fcs_join_obj,
  report_algorithm = c("leiden", "flowsom", "louvain", "phenograph", "git"),
  outdir = getwd()
)
}
\arguments{
\item{fcs_join_obj}{A list with a \verb{[[algorithm]][["abundance"]]} matrix, as produced by
FCSimple::fcs_calculate_abundance().}

\item{report_algorithm}{Character; name of the abundance matrix to export. One of
\code{"leiden"}, \code{"flowsom"}, \code{"louvain"}, \code{"phenograph"}, or \code{"git"}.}

\item{outdir}{Character; file path to an existing directory. Defaults to \code{getwd()}.}
}
\value{
Invisibly returns the abundance matrix (numeric matrix with sample rows
and cluster‐ID columns).
}
\description{
Exports the cluster‐abundance matrix (samples × clusters) to a CSV in the
specified directory, and returns the matrix invisibly.
}
\examples{
\dontrun{
  # Assume 'clustered' has abundance computed
  abundance_mat <- FCSimple::fcs_calculate_abundance(clustered)

  # Write to your working directory
  FCSimple::fcs_report_abundance(
    clustered,
    report_algorithm = "leiden",
    outdir = "~/my_analysis/results"
  )
}

}
\seealso{
FCSimple::fcs_calculate_abundance
}
