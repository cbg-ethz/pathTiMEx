% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funcsPoset.R
\name{plotPosetPaths}
\alias{plotPosetPaths}
\alias{plot_poset}
\title{Plots a poset on pathways}
\usage{
plotPosetPaths(poset, groupsInput, Datamat, size = 30)
}
\arguments{
\item{poset}{the poset to be plotted, given as a binary matrix with 
relationships between pathways}

\item{groupsInput}{the assignment of genes to pathways}

\item{Datamat}{the alterations of genes, represented as a binary matrix}

\item{size}{the size of the names of events, default to 30}
}
\value{
the plot of the input poset, together with assignment of genes to 
pathways
}
\description{
\code{plotPosetPaths} plots a poset on the level of pathways 
given as input, by representing  the causal relations between pathways,
together with gene membership.
}
\details{
The column names of the poset should be the names of the events, 
and they will be used for plotting. This function only works for groups 
generated after the stochastic search step.
}
\author{
Simona Constantinescu, \email{simona.constantinescu@bsse.ethz.ch}
}

