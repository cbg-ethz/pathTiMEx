% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcsOth.R
\name{summarizeSims}
\alias{summarizeSims}
\title{Assesses how well each group was reconstructed}
\usage{
summarizeSims(pathwaysTableSum, groups, noReps)
}
\arguments{
\item{pathwaysTableSum}{a binary matrix of all unique groups identified 
across many repetitions, together with summary information about each of the 
identified groups, as returned by the function \code{countPerfSims}.}

\item{groups}{a list with the real assignment of genes to pathways. each 
element of the list is a vector with the genes belonging to that pathway.}

\item{noReps}{the number of repetitions of the simulation experiments}
}
\value{
a list with as many elements as real pathways. each element is
a matrix in which the first column is the reconstruction percentage of
each pathway, and the second element is the number of times in which that
reconstruction was attained.
}
\description{
\code{summarizeSims} assesses the reconstruction performance
of each group, inferred from any binary matrix.
}
\details{
This function further summarizes multiple repetitions of the 
inference of pathways from a simulated dataset, having as input a summary
table with the unique groups inferred across many runs, as returned by the
function \code{countPerfSims}. The summary returned by this function doesn't
provide an assessment of the false positive rate.
}
\author{
Simona Constantinescu, \email{simona.constantinescu@bsse.ethz.ch}
}

