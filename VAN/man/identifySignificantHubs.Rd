\name{identifySignificantHubs}
\alias{identifySignificantHubs}

\title{ Identify the enriched hubs for a given combination of expression and PPI/miRnome datasets }

\description{
This function estimates the p-value for the hypothesis that the correlation between 
a hub and its interactors changes across two or more biological conditions
}

\usage{
identifySignificantHubs(exprFile, labelIndex, mapFile, outFile
		  , hubSize = 5, randomizeCount = 1000
		  , adjustMethod="BH", assocType = "TCC"
		  , labelVect = NULL
		  , exprDataType="SYMB", ppiDataType="SYMB", outputDataType="SYMB"
		  , species="Human", inputCores=4
		  , hubVect = NULL, interactomeVect = NULL)
}

\arguments{
\item{ exprFile }{ Vector of file names corresponding to normalized expression data }
\item{ labelIndex }{ Row of the exprFile which contains the sample labels }
\item{ mapFile }{ File name corresponding to PPI/Mirnome  }
\item{ outFile }{ Output file name }
\item{ hubSize }{ Minimum number of interactors in the expression dataset. Default: 5 }
\item{ randomizeCount }{ Number of permutations to consider for estimating the p-values. Default: 1000  }
\item{ adjustMethod }{ Method for adjusting the p-values. Default:"BH"}
\item{ assocType }{ Type of correlation to calculate - TCC or PCC or FSTAT. Default:TCC }
\item{ labelVect }{ Vector of conditions to test. Default: NULL }
\item{ exprDataType }{ "ENTREZ" or "SYMB". Default:SYMB }
\item{ ppiDataType }{ "ENTREZ" or "SYMB". Default:SYMB }
\item{ outputDataType }{ "ENTREZ" or "SYMB". Default:SYMB }
\item{ species }{ Name of species. At present only "Human" is supported }
\item{ inputCores }{ Number of threads for executing the code. Default:4 }
\item{ hubVect }{A vector of hub genes to consider. Default: NULL }
\item{ interactomeVect }{A vector of interactors to consider. Default: NULL }
}

\value{
Two output files are generated - one containing the p-values for the hubs and another containing
the actual correlation coefficients per condition per hub-interactor pair. While 
the name of the first output file is provided as an input parameter, the name of the second file is 
obtained by adding the suffix "Cor.txt" to the use-specified input parameter.
}

\author{ Vivek Jayaswal }

\examples{

 ## Refer User Guide
			
}

\keyword{ file }

