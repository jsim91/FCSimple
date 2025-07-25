% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_gating.R
\name{fcs_gating_object}
\alias{fcs_gating_object}
\title{Initialize Gating Object}
\usage{
fcs_gating_object(fcs_obj)
}
\arguments{
\item{fcs_obj}{A list‐like FCSimple object (e.g. from FCSimple::fcs_join) containing at
least:
\itemize{
\item data: numeric matrix or data.frame of events × channels
\item source: character vector of sample identifiers
\item run_date: (optional) acquisition dates
}}
}
\value{
An object of class “fcs_gating_object” with elements:
\itemize{
\item data
\item source
\item run_date (if present)
\item object_history: character vector recording the gating steps
}
}
\description{
Wraps a joined FCSimple object into a gating‐friendly structure.
Records the raw data matrix, sample sources, optional run dates,
and initializes an empty gating history.
}
\seealso{
FCSimple::fcs_join, FCSimple::fcs_update
}
