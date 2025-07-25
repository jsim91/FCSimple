% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fcs_project.R
\name{fcs_project_parameters}
\alias{fcs_project_parameters}
\title{Project Marker Expression onto a 2D Embedding}
\usage{
fcs_project_parameters(
  fcs_join_obj,
  override_correction = TRUE,
  reduction = c("UMAP", "tSNE"),
  parameters = "all",
  outdir = getwd(),
  sample_size = 50000,
  point_size = 0.8,
  trim_outliers = TRUE,
  trim_quantile = 0.01,
  force_xlim = FALSE,
  force_ylim = FALSE
)
}
\arguments{
\item{fcs_join_obj}{A list returned by FCSimple::fcs_join(), optionally augmented by
fcs_batch_correction() and fcs_reduce_dimensions(), containing:
\itemize{
\item \code{data}: raw or transformed expression matrix (events × channels)
\item \code{batch_correction$data} (if present)
\item a \code{umap} or \code{tsne} element with \verb{$coordinates} (events × 2)
}}

\item{override_correction}{Logical; if \code{TRUE} (default), use \code{fcs_join_obj$data} even when
batch correction exists. If \code{FALSE}, uses
\code{fcs_join_obj$batch_correction$data}.}

\item{reduction}{Character; embedding to use. Must be either \code{"UMAP"} (default) or \code{"tSNE"}.}

\item{parameters}{Character vector of channel names to project. If \code{"all"} (default),
all columns in the selected data matrix are plotted.}

\item{outdir}{Character; path to an existing directory where the PDF will be saved.
Defaults to \code{getwd()}.}

\item{sample_size}{Integer; maximum number of events to sample per parameter before plotting
(default 50000).}

\item{point_size}{Numeric; point size for the rasterized scatter (\code{geom_point_rast})
(default 0.8).}

\item{trim_outliers}{Logical; if \code{TRUE} (default), trims the lower and upper
\code{trim_quantile} fraction of observations before sampling.}

\item{trim_quantile}{Numeric; quantile threshold for outlier trimming (default 0.01).}

\item{force_xlim}{Numeric vector of length 2 to fix the x‐axis limits, or \code{FALSE}
(default) to use the data range.}

\item{force_ylim}{Numeric vector of length 2 to fix the y‐axis limits, or \code{FALSE}
(default) to use the data range.}
}
\value{
Invisibly returns \code{NULL}. The side effect is a multi‐page PDF written
to \code{outdir}.
}
\description{
Creates a multi‐panel PDF of channel (parameter) expression projected onto
a 2D reduction embedding (UMAP or tSNE). For each specified parameter,
events are optionally trimmed of outliers, subsampled, and plotted as a
rasterized scatter colored by expression intensity. Panels are arranged
in 2×2 grids and written to disk.
}
\details{
\enumerate{
\item Determines which expression matrix to use based on
\code{override_correction}.
\item Extracts UMAP or tSNE coordinates from
\code{fcs_join_obj[[ tolower(reduction) ]][["coordinates"]]}.
\item Selects the requested \code{parameters} (all by default).
\item For each parameter:
\itemize{
\item Trims outliers at the specified quantiles.
\item Subsamples up to \code{sample_size} points.
\item Creates a rasterized scatter plot colored by expression
(\code{scale_color_viridis(option="D")}).
\item Honors any fixed \code{force_xlim}/\code{force_ylim}.
}
\item Arranges plots into pages of 2×2 panels via ggpubr::ggarrange().
\item Saves a timestamped PDF named
\verb{<reduction>_parameter_projections_<YYYY-MM-DD_HHMMSS>.pdf}.
}
}
\examples{
\dontrun{
  # Assume joined, reduced and batch‐corrected object
  joined <- FCSimple::fcs_join(files)
  reduced <- FCSimple::fcs_reduce_dimensions(joined, method = "UMAP")
  corrected <- FCSimple::fcs_batch_correction(reduced)

  # Project all channels onto UMAP
  FCSimple::fcs_project_parameters(
    corrected,
    reduction = "UMAP",
    outdir    = "~/results"
  )

  # Project only CD3 and CD19
  FCSimple::fcs_project_parameters(
    corrected,
    parameters = c("CD3","CD19"),
    sample_size = 30000
  )
}

}
\seealso{
FCSimple::fcs_reduce_dimensions, FCSimple::fcs_plot_reduction,
ggplot2::ggplot, viridis::scale_color_viridis,
ggrastr::geom_point_rast, ggpubr::ggarrange
}
