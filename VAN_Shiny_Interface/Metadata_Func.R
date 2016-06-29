###################################
## Functions for handling metadata
###################################

## LICENSE:
## Copyright (C) <2012>  <Vivek Jayaswal>
## 
## This library is free software; you can redistribute it and/or modify it 
## under the terms of the GNU Lesser General Public License as published by 
## the Free Software Foundation; either version 2.1 of the License, or (at 
## your option) any later version.
## 
## This library is distributed in the hope that it will be useful, but 
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
## or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public 
## License for more details.
## 
## You should have received a copy of the GNU Lesser General Public License 
## along with this library; if not, write to the Free Software Foundation Inc., 
## 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 


#######################################################################
## Generate the mapping of gene symbol to Entrez ID using Bioconductor 
## annotation files
##
## Input
##	speciesName: Currently only "Human"
##
## Output
##	A global matrix "geneSymbEntrezMatrix" with first column 
##	corresponding to Entrez ID and second column to Gene Symbols
#######################################################################

generateGeneSymEntrezMap = function(speciesName) {
	
	geneSymbInfo = NULL
	
	if(speciesName == "Human") {
		library(org.Hs.eg.db)
		geneSymbInfo = AnnotationDbi::as.list(org.Hs.egSYMBOL)
	}
	
	## Save the Gene Symbol to Entrez ID mapping as a global variable
	geneSymbEntrezMatrix <<- matrix(names(geneSymbInfo), ncol=1)
	geneSymbEntrezMatrix <<- cbind(geneSymbEntrezMatrix, unlist(geneSymbInfo))
		
}

#######################################################################
## Obtain the Entrez Ids for a given vector of gene symbols
##
## Input
##	inputVect: Vector of gene smbols
##
## Output
##	A vector of corresponding Entrez Ids
#######################################################################

geneSymbolToEntrez = function(inputVect) {

	matchIndexes = match(inputVect, geneSymbEntrezMatrix[, 2])
	return(trimWhiteSpace(geneSymbEntrezMatrix[matchIndexes, 1]))

}

#######################################################################
## Obtain the Gene Symbols for a given vector of Entrez Ids
##
## Input
##	inputVect: Vector of Entrez Ids
##
## Output
##	A vector of corresponding gene symbols
#######################################################################

entrezToGeneSymbol = function(inputVect) {

	matchIndexes = match(inputVect, geneSymbEntrezMatrix[, 1])
	return(trimWhiteSpace(geneSymbEntrezMatrix[matchIndexes, 2]))

}

#############################################################################
## Generate the mapping of enriched hub genes to known cancer annotation
##
## Input
##	hubFile: Output file corresponsing to an expression-PPI combination
##	cancerAnnotationFile: XL-sheet containing biological annotation
##	outFile: Output file
##
## Output
##	A tab-separated file containing the CIC annotation for enriched 
##	hubs (unadjusted p-value < 0.05)
#############################################################################

obtainCancerInfo = function(hubFile, cancerAnnotationFile, outFile) {

	hubMatrix = as.matrix(read.table(hubFile, sep="\t", header=TRUE))
	lowProbIndexes = which(as.numeric(hubMatrix[, 4]) < 0.05)
	
	if(length(lowProbIndexes) == 0) return("No enriched hubs found")
	
	lowProbMatrix = hubMatrix[lowProbIndexes, ]
	lowProbHubs = lowProbMatrix[, 1]
	
	cancerMatrix = as.matrix(read.xlsx(cancerAnnotationFile, sheetIndex=1))
	cicSymb = cancerMatrix[, 1]
	
	lowProbHubsIndexes = match(lowProbHubs, cicSymb)
	foundHubs = which(!is.na(lowProbHubsIndexes))
	foundIndexes = lowProbHubsIndexes[foundHubs]
	
	numCols = ncol(cancerMatrix)
	
	if(length(foundHubs) > 0) {
	
		if(length(foundHubs) == 1) tempMatrix = matrix(cancerMatrix[foundIndexes, ], nrow=1)
		else tempMatrix = cancerMatrix[foundIndexes, ]
		
		tempMatrix = cbind(tempMatrix, lowProbMatrix[foundHubs, 4], lowProbMatrix[foundHubs, 5])
		matrixCols = c(colnames(cancerMatrix), "P-value", "Adjusted P-value")
		
		printOrder = c(1, numCols+1, numCols+2, 2:numCols)

		revMatrix = tempMatrix[, printOrder]
		if(nrow(tempMatrix) == 1) revMatrix = matrix(revMatrix, nrow=1)
		
		colnames(revMatrix) = matrixCols[printOrder]
		
		write(colnames(revMatrix), outFile, sep="\t", ncolumns=ncol(revMatrix))
		write(t(revMatrix), outFile, sep="\t", ncolumns=ncol(revMatrix), append=TRUE)
	
	}
	else print("No hub with p < 0.05 found in CIC")

}

#######################################################################
## Generate the PPI dataset from miTab lite files
##
## Input
##	mapFile: A miTab lite file
##	outFile: Output file name
##
## Output
##	Two PPI datasets - one containing the information in terms 
## 	of Entrez Ids and the second in terms of gene symbols
#######################################################################

generatePpiMap = function(mapFile, outFileHeader) {

	## Load the Entrez -> Gene Symbol mapping
	generateGeneSymEntrezMap("Human")

	## Read miTab file
	miTabMatrix = as.matrix(read.table(mapFile, header=TRUE, sep="\t"))
	
	## Obtain the Uniprot to Entrez ID mapping
	uniprotMatrix = loadFile("Uniprot", FALSE)
	
	hubIndexes = match(miTabMatrix[, 2], uniprotMatrix[, 1])
	interactorIndexes = match(miTabMatrix[, 3], uniprotMatrix[, 1])
	
	outMatrix = matrix(uniprotMatrix[hubIndexes, 3], ncol=1)
	outMatrix = cbind(outMatrix, uniprotMatrix[interactorIndexes, 3])
	
	## Identify the indexes to remove -- include NAs and blanks	
	naIndexes = which(is.na(hubIndexes))
	naIndexes = union(naIndexes, which(is.na(interactorIndexes)))
	naIndexes = union(naIndexes, which(outMatrix[, 1] == ""))
	naIndexes = union(naIndexes, which(outMatrix[, 2] == ""))			
	
	if(length(naIndexes) > 0) {
		print(paste("Number of miTab rows not found = ", length(naIndexes), sep=""))
		outMatrix = outMatrix[-naIndexes, ]
	}	
	
	## Remove semi-colons
	outMatrix = cleanData(outMatrix)
	
	## Obtain Hub-Interactor mapping based on Entrez ID
	outFileEntrez = paste(gsub(".txt", "", outFileHeader), "_Entrez.txt", sep="")
	writeHubInteractorFile(outMatrix, TRUE, outFileEntrez)
	
	## Obtain Gene Symb based mapping
	geneSymbHubs = entrezToGeneSymbol(outMatrix[, 1])
	geneSymbInteractors = entrezToGeneSymbol(outMatrix[, 2])
	
	## Identify the indexes to remove
	naIndexesSymb = which(is.na(geneSymbHubs))
	naIndexesSymb = union(naIndexesSymb, which(is.na(geneSymbInteractors)))
	
	outMatrixGene = matrix(geneSymbHubs, ncol=1)
	outMatrixGene = cbind(outMatrixGene, geneSymbInteractors)
		
	if(length(naIndexesSymb) > 0) {
		write(t(outMatrix[naIndexesSymb, ]), "Error_PPI_Generation.txt", ncolumns=2, sep="\t")	
		outMatrixGene = outMatrixGene[-naIndexesSymb, ]
	}	
	
	## Obtain Hub-Interactor mapping based on gene symbol
	outFileGeneSymb = paste(gsub(".txt", "", outFileHeader), "_Symb.txt", sep="")
	writeHubInteractorFile(outMatrixGene, TRUE, outFileGeneSymb)

}

#######################################################################
## Clean up the data in Uniprot file to ensure that a uniprot ID 
## maps to a single Entrez ID
##
## Input
##	inMatrix: A 2-column matrix with both columns containing 
##		  the Entrez IDs as defined in the Uniprot file
##
## Output
##	A 2-column output matrix with no semi-colon separate Entrez Ids
#######################################################################

cleanData = function(inMatrix) {

	tempMatrix = inMatrix
	semiColon1 = grep(";", inMatrix[, 1])
	semiColon2 = grep(";", inMatrix[, 2])
	
	changeVect1 = changeVect2 = NULL
	
	if(length(semiColon1) > 0) {		
		for(i in 1:length(semiColon1)) changeVect1 = c(changeVect1, unlist(strsplit(inMatrix[semiColon1[i], 1], split=";"))[1])
	}

	if(length(semiColon2) > 0) {		
		for(i in 1:length(semiColon2)) changeVect2 = c(changeVect2, unlist(strsplit(inMatrix[semiColon2[i], 2], split=";"))[1])
	}

	if(!is.null(changeVect1)) tempMatrix[semiColon1, 1] = changeVect1
	if(!is.null(changeVect2)) tempMatrix[semiColon2, 2] = changeVect2
	
	return(tempMatrix)

}

#######################################################################
## Generate the miRnome dataset for TargetScan or MicroCosm
##
## Input
##	databaseName: Name of miRNA-mRNA mapping database
##			Currently supports "Targetscan" and "Microcosm"
##	outFileHeader: Output file name
##
## Output
##	Two miRnome datasets - one containing the information in terms 
## 	of Entrez Ids and the second in terms of gene symbols
#######################################################################

generateMicroRnaMap = function(databaseName, outFileHeader) {

	## Load the Entrez -> Gene Symbol mapping
	generateGeneSymEntrezMap("Human")

	## Load the relavant matrix
	microMapMatrix = loadFile(databaseName, FALSE)
		
	## Obtain Hub-Interactor mapping based on Entrez ID
	outFileEntrez = paste(gsub(".txt", "", outFileHeader), "_Entrez.txt", sep="")
	writeHubInteractorFile(microMapMatrix, FALSE, outFileEntrez)
	
	## Obtain Gene Symb based mapping for the target genes
	geneSymbInteractors = entrezToGeneSymbol(microMapMatrix[, 2])
	
	## Identify the indexes to remove
	naIndexesSymb = which(is.na(geneSymbInteractors))
	
	outMatrixGene = matrix(microMapMatrix[, 1], ncol=1)
	outMatrixGene = cbind(outMatrixGene, geneSymbInteractors)
		
	if(length(naIndexesSymb) > 0) {
		write(t(microMapMatrix[naIndexesSymb, ]), "Error_Mirnome_Generation.txt", ncolumns=2, sep="\t")	
		outMatrixGene = outMatrixGene[-naIndexesSymb, ]
	}	
	
	## Obtain Hub-Interactor mapping based on gene symbol
	outFileGeneSymb = paste(gsub(".txt", "", outFileHeader), "_Symb.txt", sep="")
	writeHubInteractorFile(outMatrixGene, FALSE, outFileGeneSymb)

}

#######################################################################
## Save the hub-interactor PPI dataset in a user-defined file
##
## Input
##	inMatrix: A 2-column matrix with both columns containing 
##		  the Entrez IDs as defined in the Uniprot file
##	oFile: Output file
##
## Output
##	Tab-separated hub-interactor pairs
#######################################################################

writeHubInteractorFile = function(inMatrix, dualPair, oFile) {
	
	outMatrix = NULL
	
	if(dualPair) {
		leftVect = c(inMatrix[, 1], inMatrix[, 2])
		rightVect = c(inMatrix[, 2], inMatrix[, 1])
		outMatrix = matrix(leftVect, ncol=1)
		outMatrix = cbind(outMatrix, rightVect)
	}
	else outMatrix = inMatrix

	write(c("Hub", "Interactors"), oFile, ncolumns=2, sep="\t")
	write(t(outMatrix), oFile, ncolumns=2, sep="\t", append=TRUE)

}

#######################################################################
## Obtain the relevant metadata databse or database version
##
## Input
##	objectName:Currently supports three values - "Uniprot", 
##			"Targetscan", "Microcosm"
##	versionInfo: TRUE/FALSE
##
## Output
##	Return the R object corresponding to the database of interest
#######################################################################

loadFile = function(objectName, versionInfo) {

	data(MapMetaData)
	
	tempObject = NULL

	if(!versionInfo) {
	
		if(objectName == "Uniprot") tempObject = uniprotMatrix
		if(objectName == "Targetscan") tempObject = targetscanMatrix
		if(objectName == "Microcosm") tempObject = microcosmMatrix
	
	}	
	else {

		if(objectName == "Uniprot") tempObject = uniprotVersion
		if(objectName == "Targetscan") tempObject = targetscanVersion
		if(objectName == "Microcosm") tempObject = microcosmVersion
	
	}
			
	return(tempObject)
	
}

