\name{generateMicroRnaMap}
\alias{generateMicroRnaMap}

\title{ Obtain the miRnome mapping }

\description{
This function is used to generate an output file correponding to microRNAs as hubs 
and their potential targets as interactors.
}

\usage{
generateMicroRnaMap(databaseName, outFileHeader)
}

\arguments{
\item{ databaseName }{ Name of miRNA-mRNA mapping database. Currently supports "Targetscan" and "Microcosm" }
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

