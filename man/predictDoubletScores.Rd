% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doublets.r
\name{predictDoubletScores}
\alias{predictDoubletScores}
\title{predictDoubletScores}
\usage{
predictDoubletScores(dbobj, model, reduction)
}
\arguments{
\item{dbobj}{An amethyst doublet object containing the data with artificial doublets; output of makeDoubletObject}

\item{model}{A trained model used to predict doublet probabilities; output of buildDoubletMode}

\item{reduction}{A character string specifying the name of the reduction slot to use.}
}
\value{
The amethyst doublet object with updated metadata including doublet scores.
}
\description{
This function calculates doublet scores for cells and adds them to the metadata of the provided doublet object.
}
\examples{
\dontrun{
  dbobj <- predictDoubletScores(dbobj, reduction = "irlba", model = result$model)
}
}
