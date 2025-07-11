% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dim_reduction.R
\name{dimEstimate}
\alias{dimEstimate}
\title{dimEstimate}
\usage{
dimEstimate(obj, genomeMatrices, dims, threshold = 0.98)
}
\arguments{
\item{obj}{Amethyst object containing the matrix to calculate, which should be in the genomeMatrices slot}

\item{genomeMatrices}{Name of the matrix in the genomeMatrices slot to calculate nv}

\item{dims}{Number of singular values to test for each matrix}

\item{threshold}{Amount of variance that must be explained by the nv value, with 1 being 100\% of variance}
}
\value{
Integer indicating the number of principal components required to meet the variance threshold for each matrix
}
\description{
Estimate the nv value needed for singular value decomposition with irlba
}
\examples{
\dontrun{
  dimEstimate(obj = combined, genomeMatrices = c("cg_100k_score", "ch_100k_pct"), dims = c(50, 50), threshold = 0.98)
}

}
