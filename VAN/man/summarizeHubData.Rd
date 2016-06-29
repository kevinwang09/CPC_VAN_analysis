\name{summarizeHubData}
\alias{summarizeHubData}

\title{ Combine the results obtained using multiple expression and/or PPI datasets }

\description{
This function generates an output file that summarizes the p-values obtained using 
multiple expression datasets and/or PPI datasets. The p-values are summarized using 
Fisher's combined test or Rank Prod method.
}

\usage{
summarizeHubData(fileNames=NULL, outFile, metaAnalysis=NULL, rankProdItr=NULL)
}

\arguments{
\item{ fileNames }{ Vector of file names corresponding to individual expression-PPI datasets results. If NULL, all files with the 
		    extension txt are considered. Default: NULL }
\item{ outFile }{ File name that contains the combined information }
\item{ metaAnalysis }{ Type of meta-analysis to perform
			NULL: No meta-analysis to perform.
			Fisher: Fisher's combined analysis
			RankProd: Rank-prod analysis 
			Default: NULL}
\item{ rankProdItr }{ Number of iteration to perform for significance of rank prod values. Default: NULL. 
			The default value should be used only if rank prod analysis is not being performed }
}

\value{
Summarized Output file with optional columns corresponding to meta-analysis results
}

\author{ Vivek Jayaswal }

\examples{

## Refer user guide
			
}

\keyword{ file }

