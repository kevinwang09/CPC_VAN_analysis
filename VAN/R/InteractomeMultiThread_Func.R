##################################################################
## Functions for analysing the interactome and estimating the 
## p-values in a multi-threaded setting
##################################################################

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
## For a given expression dataset and PPI dataset, estimate the p-value 
## for each hub with >=5 interaction partners
##
## Input
##	exprFile: Vector of file names corresponding to normalized expression data
##	labelIndex: Row of the exprFile which contains the sample labels
##	mapFile: File name corresponding to PPI/Mirnome
##	outFile: Output file name
##	hubSize: Minimum number of interactors in the expression dataset
##	randomizeCount: Number of permutations to consider for estimating 
##			the p-values
##	adjustMethod: Method for adjusting the p-values. Default:"BH"
##			Possible values - "BH", "bonferroni"
##	assocType: Type of correlation to calculate. Default:TCC
##		   TCC, PCC, FSTAT
##	labelVect: Vector of conditions to test. If all conditions 
##		   are to be tested, set to NULL. Default: NULL
##	exprDataType: "ENTREZ" or "SYMB". Default:SYMB
##	ppiDataType: "ENTREZ" or "SYMB". Default:SYMB
##	outputDataType: "ENTREZ" or "SYMB". Default:SYMB
##	species: Name of species. At present only "Human" is supported
##	inputCores: Number of threads for executing the code. Default:4
##	hubVect: A vector of hub genes to consider. Default: NULL
##	interactomeVect: A vector of interactors to consider. Default: NULL
##
##
## Output
##	Two output files - 
##	1. outFile that contains the p-value for each hub
##	2. outFile with suffix "Cor.txt" that contains the CC value
##	   for each hub-interactor pair
#######################################################################

identifySignificantHubs = function(exprFile, labelIndex, mapFile, outFile
				  , hubSize = 5, randomizeCount = 1000
				  , adjustMethod="BH", assocType = "TCC"
				  , labelVect = NULL
				  , exprDataType="SYMB", ppiDataType="SYMB", outputDataType="SYMB"
				  , species="Human", inputCores=4
				  , hubVect = NULL, interactomeVect = NULL) {
	
	## Create a multi-threaded system
	numCores = inputCores
	if(detectCores() < inputCores) numCores = detectCores()

	cl = makeCluster(numCores)
	registerDoParallel(cl)
	
	twoInputExprFiles = FALSE
	if(length(exprFile) == 2) twoInputExprFiles = TRUE
	
	## Load the Bioconductor annotation file
	generateGeneSymEntrezMap(species)
	
	## Read expression data for genes and check the association type
	srcExprMatrix = readExprData(exprFile[1], labelIndex)
	
	## Ensure that the correct statistic is being calculated 
	## and all the labels are present in the expression matrix
	labelsToConsider = colnames(srcExprMatrix)
	if(!is.null(labelVect)) labelsToConsider = labelVect
	
	uniqueLabels = unique(labelsToConsider)
	correctAssocType = checkAssociation(uniqueLabels, assocType)
	stopifnot(correctAssocType == TRUE)
		
	## Read mapping data
	srcHubsMatrix = trimWhiteSpace(as.matrix(read.table(mapFile, sep="\t", header=TRUE)))
	
	## Perform an internal conversion to Entrez IDs if expr and ppi data types 
	## are different
	if(exprDataType != ppiDataType) {
	
		naIndexesExpr = naIndexesPpiHubs = naIndexesPpiInt = NULL
	
		if(exprDataType != "ENTREZ") {

			entrezIdVect = geneSymbolToEntrez(rownames(srcExprMatrix))						
			naIndexesExpr = which(is.na(entrezIdVect))
			naGeneSymbols = rownames(srcExprMatrix)[naIndexesExpr]
			rownames(srcExprMatrix) = entrezIdVect
			
			if(length(naIndexesExpr) > 0) {			
				print("Missing Entrez IDs for expression data")
				generateErrorOutput(naGeneSymbols, "Expr")
				srcExprMatrix = srcExprMatrix[-naIndexesExpr, ]
			}

		}
		
		if(ppiDataType != "ENTREZ") {

			if(!twoInputExprFiles) { 
			
				## If two files are provided, then hubs are microRNAs and no conversion needed for hubs				
				entrezIdVect = geneSymbolToEntrez(srcHubsMatrix[, 1])					
				naIndexesPpiHubs = which(is.na(entrezIdVect))
				naGeneSymbols = srcHubsMatrix[naIndexesPpiHubs, 1]
				srcHubsMatrix[, 1] = entrezIdVect
				
				if(length(naIndexesPpiHubs) > 0) {				
					print("Missing Entrez IDs in hubs")
					generateErrorOutput(naGeneSymbols, "PPI_Hubs")
				}
			
			}	
			
			entrezIdVect = geneSymbolToEntrez(srcHubsMatrix[, 2])
			naIndexesPpiInt = which(is.na(entrezIdVect))
			naGeneSymbols = srcHubsMatrix[naIndexesPpiInt, 2]
			srcHubsMatrix[, 2] = entrezIdVect
			
			if(length(naIndexesPpiInt) > 0) {				
				print("Missing Entrez IDs in interactors")
				generateErrorOutput(naGeneSymbols, "PPI_Interactors")
			}
			
		}
		
		if(!is.null(naIndexesPpiHubs) || !is.null(naIndexesPpiInt)) srcHubsMatrix = srcHubsMatrix[-c(naIndexesPpiHubs, naIndexesPpiInt), ]
		
		if(is.null(naIndexesExpr) || is.null(naIndexesPpiHubs) || is.null(naIndexesPpiInt)) {
			print("Do you want to continue: (Y/N)? ")

			userInput = "X"			
			while(is.na(match(userInput, c("Y", "N")))) userInput = readLines(n=1)
			
			if(userInput == "N") return(0)
		}	
				
	}

	## Obtain the list of hubs
	srcHubs = readMapData(srcHubsMatrix)
	elementsInHubs = sapply(srcHubs, length)

	## Read regulator expression data
	srcRegMatrix = srcExprMatrix
	if(twoInputExprFiles) srcRegMatrix = readExprData(exprFile[2], labelIndex)	
	regNames = rownames(srcRegMatrix)

	## Identify the hubs based on expression data
	exprHubs = filterHubs(srcHubs, rownames(srcExprMatrix), regNames, hubSize)
	elementsInExprHub = sapply(exprHubs, length)

	## Obtain user-defined subset of hubs and interactomes
	if(!is.null(hubVect) || !is.null(interactomeVect)) {
		exprHubs = filterUserDefined(exprHubs, hubVect, interactomeVect, hubSize)
		elementsInExprHub = sapply(exprHubs, length)
	}	
	
	exprIndexes = match(names(exprHubs), names(srcHubs))
	totalElementsInExprHub = elementsInHubs[exprIndexes]

	## Save the hubs/interactors as gene symbols if the user explicitly asks for it
	## and the input comprises Entrez IDs
	changeToGeneSymb = FALSE
	if(outputDataType == "SYMB" && (exprDataType == "ENTREZ" || ppiDataType == "ENTREZ")) changeToGeneSymb = TRUE

	## Display the number of threads used for perumation tests
	print(paste("Number of threads = ", getDoParWorkers(), sep=""))
	
	## Perform interactome analysis
	filterSrcExprMatrix = filterMatrix(srcExprMatrix, uniqueLabels)
	filterSrcRegMatrix = filterMatrix(srcRegMatrix, uniqueLabels)
	
	probVect = interactomeAnalysis(filterSrcExprMatrix, filterSrcRegMatrix, exprHubs, randomizeCount, assocType, numCores
					, outFile, changeToGeneSymb, twoInputExprFiles, uniqueLabels)	
					
	saveProbValues(probVect, elementsInExprHub, totalElementsInExprHub, outFile, changeToGeneSymb, twoInputExprFiles, adjustMethod)
	
	stopCluster(cl)

}

#######################################################################
## For a given normalized expression matrix and PPI list, 
## estimate the p-value for each hub with >=5 interaction partners
##
## Input
##	inExprMatrix: Normalized expression matrix for interactors
##	inRegMatrix: Normalized expression matrix for regulators
##	inHubList: PPI list with elements corresponding to the interactors
##	inCount: Number of permutations to consider for estimating 
##		 the p-values
##	assocType: Type of correlation to calculate
##	coreCount: Number of threads for executing the code
##	fileName: Output file name
##	getGeneSymb: TRUE/FALSE. If TRUE, convert Entrez IDs to gene symbols
##	isMicro: TRUE/FALSE. If TRUE, then the regulator corresponds to microRNAs
##	sampleGroups: Vector of unique conditions to evaluate
##
##
## Output
##	A vector of p-values for each hub
#######################################################################

interactomeAnalysis = function(inExprMatrix, inRegMatrix, inHubList, inCount, assocType, coreCount, fileName
				, getGeneSymb, isMicro, sampleGroups) {

	probVect = NULL

	allLabels = colnames(inExprMatrix)
	numUniqueLabels = length(sampleGroups)

	hubNames = names(inHubList)
	geneNames = rownames(inExprMatrix)
	regulatorNames = rownames(inRegMatrix)

	## Create the file for saving <hub, interactor> values for the two conditions	
	corFile = gsub(".txt", "_Cor.txt", fileName)
	
	if(numUniqueLabels == 2) {
		corLabel = paste(sampleGroups[1], sampleGroups[2], sep="-")
		write(c("Hub", "Interactor", sampleGroups, corLabel), corFile, sep="\t", ncolumns=5)
	}
	else write(c("Hub", "Interactor", sampleGroups), corFile, sep="\t", ncolumns=numUniqueLabels+2)

	print(paste("Number of hubs = ", length(inHubList), sep=""))
	
	for(hubIndex in 1:length(inHubList)) {
		
		print(paste("Current hub = ", hubIndex, sep=""))
		
		currentHub = hubNames[hubIndex]
		currentInteractors = as.character(inHubList[[hubIndex]])
		numInteractors = length(currentInteractors)
		
		currentHubIndex = match(currentHub, regulatorNames)
		currentIteratorIndexes = match(currentInteractors, geneNames)
					
		## Obtain the avg correlation coefficient for real data
		## along with the correlation coefficient per interactor
		origHubDiff = getCorrelation(inRegMatrix[currentHubIndex, ], inExprMatrix[currentIteratorIndexes, ]
					     , sampleGroups, corType = assocType, corInfo = TRUE, permuteLabels = FALSE)
		
		saveCorValues(currentHub, currentInteractors, origHubDiff$corMatrix
				, corFile, getGeneSymb, isMicro)
		
		## Obtain the correlation coefficient for randomized data
		grpValue = floor(inCount/coreCount)
		grpValueVector = rep(grpValue, coreCount)
		if(sum(grpValueVector) < inCount) grpValueVector[coreCount] = grpValueVector[coreCount] + (inCount - sum(grpValueVector))
		
		exportFunctions = c("calculateProbDist", "getCorrelation"
				, "interactomeTaylorCorrelation", "getTaylorCor"
				, "interactomePearsonCorrelation", "getPearsonCor"
				, "getSampleIndexes", "getFStat")
					
		probTemp = foreach(i = 1:coreCount, .combine='c', .export=exportFunctions) %dopar% {
		
			calculateProbDist(inExprMatrix, inRegMatrix, currentHubIndex, currentIteratorIndexes
						, sampleGroups, assocType
						, grpValueVector[i], origHubDiff$avgData)
												
		}		

		probVect = c(probVect, sum(probTemp)/inCount)		
	
	} ## All hubs have been considered
	
	names(probVect) = names(inHubList)
	return(probVect)

}

#######################################################################
## For a given hub and its interactor set, estimate the p-values
##
## Input
##	exprMatrix: Normalized expression matrix for interactors
##	regMatrix: Normalized expression matrix for regulators
##	hubIndex: Row index of regMatrix
##	interactorIndexes: Vector of row indexes in exprMatrix corresponding 
##			to the interactors
##	sampleGroups: Vector of unique conditions to compare
##	assocType: Type of correlation to calculate
##	grpVal: Number of permutations
##	origVal: Actual avg hub diff
##
##
## Output
##	Number of times the bootstrap p-value > origVal
#######################################################################

calculateProbDist = function(exprMatrix, regMatrix, hubIndex, interactorIndexes, sampleGroups
			, assocType, grpVal, origVal) {

	randomHubDiff = 0	

	for(j in 1:grpVal) {

		temp = getCorrelation(regMatrix[hubIndex, ], exprMatrix[interactorIndexes, ], sampleGroups
				, corType = assocType, corInfo = FALSE, permuteLabels = TRUE)
				
		if(temp >= origVal) randomHubDiff = randomHubDiff + 1

	}

	return(randomHubDiff)
	
}		

#######################################################################
## For a given hub, calculate the correlation for each hub-interactor pair 
## and the average hub difference
##
## Input
##	inVect: Vector of expression value for the hub
##	inMatrix: Matrix of expression values for the interactors
##	sampleGroups: Vector of unique conditions to compare
##	corType: Type of correlation to calculate
##	corInfo: TRUE -> Return a list of values
##		 FALSE -> Return only the average hub difference
##	permuteLabels: TRUE/FALSE. If TRUE, then permute the samples for 
##			calculation of p-value
##
## Output
##	If corInfo = TRUE, a list containing average hub difference and the 
##	pairwise hub-interactor values for both the conditions
##	If corInfo = FALSE, average hub difference
#######################################################################

getCorrelation = function(inVect, inMatrix, sampleGroups, corType, corInfo, permuteLabels) {

	numInteractors = nrow(inMatrix)
	avgDiffLabels = NULL
	
	corVectMatrix = matrix(0, nrow=nrow(inMatrix), ncol=length(sampleGroups))	
	labelList = getSampleIndexes(sampleGroups, colnames(inMatrix), permuteLabels)

	if(corType == "TCC") { ## Taylor's CC		
		for(i in 1:2) corVectMatrix[, i] = interactomeTaylorCorrelation(inVect, inMatrix, labelList[[i]])
	}	
	
	if(corType == "PCC" || corType == "FSTAT") { ## Pearsons CC
		for(i in 1:length(sampleGroups)) corVectMatrix[, i] = interactomePearsonCorrelation(inVect[labelList[[i]]], inMatrix[, labelList[[i]]])
	}	

	if(corType == "TCC" || corType == "PCC") {
	
		avgDiffLabels = sum(abs(corVectMatrix[, 1] - corVectMatrix[, 2]))
		avgDiffLabels = avgDiffLabels/(numInteractors - 1)
	
	}	

	if(corType == "FSTAT") avgDiffLabels = getFStat(corVectMatrix)

	if(corInfo) return(list(avgData = avgDiffLabels, corMatrix = corVectMatrix))

	return(avgDiffLabels)
}

#######################################################################
## Obtain the Taylor's CC for all interactors for a given condition
##
## Input
##	hubVector: Vector of expression value for the hub
##	geneMatrix: Matrix of expression values for the interactors
##	currentLabel: Column indexes corresponding to a condition
##
## Output
##	A vector of CC for all hub-interactor pairs
#######################################################################

interactomeTaylorCorrelation = function(hubVector, geneMatrix, currentLabel) {		
	return(apply(geneMatrix, 1, getTaylorCor, hubVector, currentLabel))
}

#######################################################################
## Obtain the Taylor's CC for a given hub-interactor pair
##
## Input
##	x: Vector of expression value for the interactor
##	hubVector: Vector of expression value for the hub
##	currentLabels: Column indexes corresponding to a condition
##
## Output
##	Taylor's CC
#######################################################################

getTaylorCor = function(x, hubVector, currentLabels) {

	numSamples = length(currentLabels)
	sdHub = sd(hubVector[currentLabels])
	sdInteractor = sd(x[currentLabels])
	denom = (numSamples - 1) * sdHub * sdInteractor

	numValue1 = x[currentLabels] - mean(x)
	numValue2 = hubVector[currentLabels] - mean(hubVector)
	ratioVal = sum(numValue1 * numValue2)/denom
	
	return(ratioVal)
}

#######################################################################
## Obtain the Pearsons CC for all interactors for a given condition
##
## Input
##	hubVector: Vector of expression value for the hub
##	geneMatrix: Matrix of expression values for the interactors
##
## Output
##	A vector of CC for all hub-interactor pairs
#######################################################################

interactomePearsonCorrelation = function(hubVector, geneMatrix) {	
	return(apply(geneMatrix, 1, getPearsonCor, hubVector))	
}

#######################################################################
## Obtain the Pearsons CC for a given hub-interactor pair
##
## Input
##	x: Vector of expression value for the interactor
##	inVect: Vector of expression value for the hub
##
## Output
##	Pearsons CC
#######################################################################

getPearsonCor = function(x, inVect) {
	return(cor(x, inVect))
}

#######################################################################
## Check that the association measure is appropriate for analyzing the 
## number of conditions in the dataset
##
## Input
##	uniqueLabels: Vector of unique conditions in the dataset
##	assocType: Type of association
##
## Output
## 	TRUE/FALSE value. If the assocType is "TCC" or "PCC", then 
## 	the number of unique conditions has to be two. Otherwise, the 
## 	the number of unique conditions must be > 2
#######################################################################

checkAssociation = function(uniqueLabels, assocType) {
	
	numUniqueLabels = length(uniqueLabels)

	stopifnot(numUniqueLabels > 1)

	if(numUniqueLabels == 2 && is.null(match(assocType, c("TCC", "PCC")))) {
		
		print("Number of unique labels = 2. The valid options are TCC/PCC")
		return(FALSE)
	}
	
	if(numUniqueLabels > 2 && assocType != "FSTAT") {
	
		print("Number of unique labels > 2. The only valid option is FSTAT")
		return(FALSE)
	}	
	
	return(TRUE)

}

#######################################################################
## Obtain a subset of the expression matrix corresponding to the 
## conditions of interest
##
## Input
##	inputMatrix: Input matrix (N x P)
##	labelsOfInterest: Vector of biological conditions of interest
##
## Output
##	An N x P1 matrix such that only the columns corresponding to
##	the relevant conditions are retained
#######################################################################

filterMatrix = function(inputMatrix, labelsOfInterest) {

	allLabels = colnames(inputMatrix)
	colsOfInterest = which(allLabels %in% labelsOfInterest)
	revMatrix = inputMatrix[, colsOfInterest]
	
	return(exprDataStd(revMatrix))

}

#######################################################################
## Convert a gene expression matrix with multiple rows corresponding
## to the same gene into a normalized matrix with one row per gene
## Also, the gene expression values are standardized with row median
## set to 0 and row var = 1
##
## Input
##	summaryMatrix: Non-standardized input matrix
##
## Output
##	Standardized matrix
#######################################################################

exprDataStd = function(summaryMatrix) {

	## Median center the values and set variance to 1
	medianVect = apply(summaryMatrix, 1, median)
	sdVect = apply(summaryMatrix, 1, sd)
	
	normalizedMatrix = summaryMatrix - medianVect
	normalizedMatrix = normalizedMatrix/sdVect

	return(normalizedMatrix)
}

#######################################################################
## Determine the column indexes of samples that correspond to different 
## biological conditions of interest
##
## Input
##	uniqueLabels: Vector of unique conditions
##	allLabels: Original labels for the various columns
##	permuteLabels: TRUE/FALSE
##
## Output
##	A list with each element representing the samples that correspond
##	to the biological condition of interest. The order of list 
##	elements corresponds to uniqueLabels
#######################################################################

getSampleIndexes = function(uniqueLabels, allLabels, permuteLabels) {

	numUniqueLabels = length(uniqueLabels)	
	
	labelIndexesList = list()	
	for(i in 1:numUniqueLabels) labelIndexesList[[i]] = which(allLabels %in% uniqueLabels[i])
	
	if(!permuteLabels) return(labelIndexesList)
	
	## Proceed if labels to be permuted
	countPerGrp = sapply(labelIndexesList, length)
	totalCount = sum(countPerGrp)
	
	grpData = NULL
	for(i in 1:numUniqueLabels) grpData = c(grpData, rep(i, countPerGrp[i]))
	
	temp = sample(grpData, totalCount, replace=FALSE)	
	labelIndexesList = list()
	for(i in 1:numUniqueLabels) labelIndexesList[[i]] = which(temp == i)
	
	return(labelIndexesList)
	
}

#######################################################################
## Determine the ratio of between to within sum of squares for testing 
## whether there is a difference between three or more conditions
## for a given hub
##
## Input
##	inputMatrix: X x Y matrix where X corresponds to the TCC/PCC
##		     value for hub-interactor pairs and Y denotes 
##		     the number of conditions
##
## Output
##	Ratio 
#######################################################################

getFStat = function(inputMatrix) {
	 
	inputData = as.vector(inputMatrix)
	totalSumSquare = var(inputData) * (length(inputData) - 1)
	
	betweenSumSquare = sum(apply(inputMatrix, 2, var) * (nrow(inputMatrix) - 1))
	withinSumSquare = totalSumSquare - betweenSumSquare 
	ratio = betweenSumSquare/withinSumSquare
	
	return(ratio)
	
}

