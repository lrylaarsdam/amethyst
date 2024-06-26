% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{makeFuzzyGeneWindows}
\alias{makeFuzzyGeneWindows}
\title{makeFuzzyGeneWindows}
\usage{
makeFuzzyGeneWindows(obj, type, ref, threads = 1)
}
\arguments{
\item{obj}{Amethyst object for which to calculate \% methylation over all genes}

\item{type}{Type of methylation to calculate; e.g. "CH".
"CG" is an option but may be less biologically relevant}

\item{ref}{Species reference to obtain bed file in 500 bp genomic windows. "hg38", "mm10", or path to bed file.}

\item{threads}{Enables multithreading}
}
\value{
Returns a matrix just like the output of makeWindows with rows as genes, columns as cells, and values as \% methylation
}
\description{
Rapidly calculate \% methylation over all genes by averaging 500 bp window values.
Slightly less accurate than the alternative indexing approach for short genes.
}
\examples{
gene_ch <- makeFuzzyGeneWindows(obj = obj, type = "CH", threads = 20)
}
