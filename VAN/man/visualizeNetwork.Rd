\name{visualizeNetwork}
\alias{visualizeNetwork}

\title{ Visualize the hub-interactor pairs }

\description{
This function enables the visualization of changes in correlation between a 
hub and its interactors across two biologial conditions. 
}

\usage{
visualizeNetwork(inputFile, inputHub, paletteVector, numThresh=NULL, condVector=NULL)
}

\arguments{
\item{ inputFile }{ File containing the hub-interactor pairs and their correlation values in each condition }
\item{ inputHub }{ Name of the hub }
\item{ paletteVector }{ A vector of three colors }
\item{ numThresh }{ Number of interactors to plot. This value is set to NULL if all the interctors should be plotted. Default: NULL }
\item{ condVector }{ Vector of the two states to plot. Default is NULL }
}

\value{
A graph with two columns such that each column represent a condition of interest. Within a 
column, a user-defined color-coding scheme is used to visualize the strength of relationship 
between a hub and its interactors.

If the inputFile corresponds to an analysis of three or more conditions, then the parameter condVector should be 
specified as a vector comprising the two states of interest.
}

\author{ Vivek Jayaswal }

\examples{

## Refer User Guide

}

\keyword{ file }

