% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{addCellInfo}
\alias{addCellInfo}
\title{addCellInfo}
\usage{
addCellInfo(obj, file)
}
\arguments{
\item{obj}{Amethyst object to add cell info metadata}

\item{file}{Path of the cellInfo.txt file}
}
\value{
Returns the same amethyst object with information contained in the cellInfo file added to the obj@metadata slot
}
\description{
Add the information contained in the cellInfo file outputs to the obj@metadata slot
}
\examples{
obj <- addCellInfo(obj = obj, file = "~/Downloads/cellInfo.txt")
}
