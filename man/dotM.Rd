% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{dotM}
\alias{dotM}
\title{dotM}
\usage{
dotM(
  obj,
  sizeMatrix,
  colorMatrix,
  genes,
  groupBy,
  splitBy = NULL,
  nrow = NULL,
  colors = NULL
)
}
\arguments{
\item{obj}{The amethyst object to plot}

\item{sizeMatrix}{Name of matrix to plot values as "size" aesthetic}

\item{colorMatrix}{Name of matrix to plot values as "color" aesthetic}

\item{genes}{List of genes to plot methylation on the x axis}

\item{groupBy}{A categorical variable contained in the metadata to group cells by on the y axis}

\item{splitBy}{Optional additional facet level}

\item{nrow}{If splitBy is specified, this option controls the number of rows the facets are distributed across}

\item{colors}{Optional; specify color palette (default is viridis)}
}
\value{
Returns a ggplot object displaying average methylation values for each gene by the grouping variable
}
\description{
Generates a dot plot of percent methylation over a gene body
}
\examples{
dotM( obj, sizeMatrix = "gene_ch", colorMatrix = "gene_ch_norm", groupBy = "type", splitBy = "batch", genes = c("SATB2", "ADARB2", "GAD1"))
}
