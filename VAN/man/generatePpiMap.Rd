\name{generatePpiMap}
\alias{generatePpiMap}

\title{ Obtain the PPI mapping }

\description{
This function is used to generate an output file correponding to protein-protein interaction. 
The input should be a MINI-Tab file downlaoded from Wodak lab website
}

\usage{
generatePpiMap(mapFile, outFileHeader)
}

\arguments{
\item{ mapFile }{ MINI tab file }
\item{ outFileHeader }{ Output file name }
}

\value{
Two tab-separated output files. The first file saves the target genes (or interactors) as 
Entrez IDs and the second saves the interactors as gene symbols. The first file has the suffix 
"Entrez.txt" and the latter has the suffix "Symb.txt"
}

\author{ Vivek Jayaswal }

\examples{

## Refer User Guide

}

\keyword{ file }

