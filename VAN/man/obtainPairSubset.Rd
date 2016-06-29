\name{obtainPairSubset}
\alias{obtainPairSubset}

\title{ Obtain a subset of the hub-interactor pairs file }

\description{
This function enables the extraction of only those hub-interactor pairs
that correspond to (a) a pre-determined set of hubs or (b) significant hubs, 
i.e. hubs with p-values below a threshold
}

\usage{
obtainPairSubset(filePrefix, useAdjustedProb, probThresh, hubNames=NULL)
}

\arguments{
\item{ filePrefix }{ File name prefix corresponding to the enriched hubs file and hub-interactor pairs file }
\item{ useAdjustedProb }{ Boolean value. TRUE -> Consider adjusted p-values. FALSE -> Consider unadjusted p-values }
\item{ probThresh }{ p-value cut-off }
\item{ hubNames }{ Vector of hub names. Default: NULL }
}

\value{
An output file with the name <filePrefix>_Cor_Signif.txt.

}

\author{ Vivek Jayaswal }

\examples{

## Refer User Guide

}

\keyword{ file }

