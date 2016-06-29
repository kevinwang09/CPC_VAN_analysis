\name{obtainCancerInfo}
\alias{obtainCancerInfo}

\title{ Identify the hubs that map to known cancer genes }

\description{
This function enables the identification of genes that act as hubs 
and have previously been shown to be associated with cancer. The latter 
information is downloaded as an XL-sheet from the Sanger institute's website.
}

\usage{
obtainCancerInfo(hubFile, cancerAnnotationFile, outFile)
}

\arguments{
\item{ hubFile }{ Output file corresponsing to an expression-PPI combination }
\item{ cancerAnnotationFile }{ XL-sheet containing biological annotation }
\item{ outFile }{ Output file }
}

\value{
A tab-separated file containing the CIC annotation for enriched hubs, i.e. hubs with unadjusted p-value < 0.05
}

\author{ Vivek Jayaswal }

\examples{

## Refer User Guide
			
}

\keyword{ file }

