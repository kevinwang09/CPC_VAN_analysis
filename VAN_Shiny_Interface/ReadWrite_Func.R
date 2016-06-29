#########################################################
## Functions for reading and writing the data
#########################################################

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
## Convert the expression dataset into a matrix
##
## Input
##	inFile: Expression dataset file
##	labelIndex: Row in inFile corresponding to sample labels
##
## Output
##	A N1 x P matrix with N1 corresponding to the number of unique genes
#######################################################################

readExprData = function(inFile, labelIndex) {
	
	tempVect = scan(inFile, sep="@", what="character")
	labelEnd = match("LABEL_END", trimWhiteSpace(gsub("\t", "", tempVect))) - 1

	## Obtain labels
	labelList = list()
	for(i in 1:labelEnd) labelList[[i]] = trimWhiteSpace(scan(inFile, what="character", nlines=1, sep="\t", skip=i-1))[-1]
	
	## Obtain gene expression data
	geneExprData = read.table(inFile, sep="\t", header=FALSE, skip=labelEnd+1)
	geneExprMatrix = as.matrix(geneExprData[, -1])	## Remove the column corresponding to Gene IDs
	rownames(geneExprMatrix) = trimWhiteSpace(geneExprData[, 1])
	colnames(geneExprMatrix) = labelList[[labelIndex]]

	return(geneExprMatrix)
}

#######################################################################
## Obtain a list with each element representing the interactors per 
## hub
##
## Input
##	mapMatrix: A matrix with the first column containing the 
##		   hubs and the second column containing the interactors
##
## Output
##	A list
#######################################################################

readMapData = function(mapMatrix) {

	hubConnectors = split(mapMatrix[, 2], mapMatrix[, 1])
	return(hubConnectors)
}

#######################################################################
## Filter the hub list so that only those hubs that contain >=5 
## interactors in the expression datset are retained
##
## Input
##	inHubList: List comprising hubs and their interactors
##	inGeneNames: Vector of genes in expression dataset
##	inRegNames: Vector of hubs in regulator dataset
##	thresh: Only consider hubs with number of interactors >= thresh
##
## Output
##	Modified list
#######################################################################

filterHubs = function(inHubList, inGeneNames, inRegNames, thresh) {

	## Hubs that are found in the regulator expression data
	exprHubIndexes = which(names(inHubList) %in% inRegNames)
	exprOnlyHubs = inHubList[exprHubIndexes]

	revHubs = list()
	revHubsIndex = 1
	revHubsNames = NULL
	currentHubNames = names(exprOnlyHubs)
	
	for(i in 1:length(exprOnlyHubs)) {
		temp = intersect(exprOnlyHubs[[i]], inGeneNames)		
		
		if(length(temp) >= thresh) {
			revHubs[[revHubsIndex]] = temp
			revHubsIndex = revHubsIndex + 1
			revHubsNames = c(revHubsNames, currentHubNames[i])
		}
		
	} ## All hubs have been evaluated
	
	names(revHubs) = revHubsNames
	return(revHubs)
}

#######################################################################
## Filter the hub list so that only those hubs and interactomes that
## are user-defined are considered
##
## Input
##	inExprHubs: List comprising hubs and their interactors
##	inHubVect: Vector of user-defined hub genes
##	inInteractomeVect: Vector of user-defined interactome genes
##	thresh: Only consider hubs with number of interactors >= thresh
##
## Output
##	Modified list
#######################################################################

filterUserDefined = function(inExprHubs, inHubVect, inInteractomeVect, thresh) {

	hubNames = names(inExprHubs)
	
	tempExprHubs = inExprHubs
	if(!is.null(inHubVect)) tempExprHubs = tempExprHubs[which(hubNames %in% inHubVect)]
	
	stopifnot(length(tempExprHubs) > 0)

	if(is.null(inInteractomeVect)) return(tempExprHubs)

	## Filter the interactomes and only consider the user-defined 
	## subset of genes as the interactome
	
	revHubs = list()
	revHubsIndex = 1
	revHubsNames = NULL
	hubNames = names(tempExprHubs)
	
	for(i in 1:length(tempExprHubs)) {
		temp = intersect(tempExprHubs[[i]], inInteractomeVect)		
		
		if(length(temp) >= thresh) {
			revHubs[[revHubsIndex]] = temp
			revHubsIndex = revHubsIndex + 1
			revHubsNames = c(revHubsNames, hubNames[i])
		}
		
	} ## All hubs have been evaluated
	
	names(revHubs) = revHubsNames
	return(revHubs)

}

#######################################################################
## Save the p-values for each hub
##
## Input
##	pVect: Vector of p-value per hub
##	elementsExpr: Number of interactors per hub in the PPI database
##			and expression dataset
##	elementsTotal: Number of interactors per hub in only the PPI database
##	fileName: Output file name
##	getGeneSymb: TRUE - Assumes that the analysis was based on Entrez 
##			    Ids and the gene symbols have to be saved explicitly 
##		     FALSE - Do nothing
##	isMicro: TRUE/FALSE. If TRUE, then the regulator corresponds to microRNAs
##	adjustMethod: Adjustment for multiple comparison
##
## Output
##	Output file
#######################################################################

saveProbValues = function(pVect, elementsExpr, elementsTotal, fileName, getGeneSymb, isMicro, adjustMethod) {

	hubNames = names(pVect)
	if(getGeneSymb && !isMicro) hubNames = entrezToGeneSymbol(hubNames)

	adjustedProb = p.adjust(pVect, method=adjustMethod)

	outMatrix = matrix(hubNames, ncol=1)
	outMatrix = cbind(outMatrix, elementsTotal)
	outMatrix = cbind(outMatrix, elementsExpr)
	outMatrix = cbind(outMatrix, pVect)
	outMatrix = cbind(outMatrix, adjustedProb)
	
	################ Original Code ####################	
# 	write(c("Hubs", "Num Interactors", "Interactors Found", "P-value", "Adjusted P-value"), fileName, ncolumns=5, sep="\t")
# 	write(t(outMatrix), fileName, ncolumns=5, sep="\t", append=TRUE)
	##################################################
	outDataFrame = as.data.frame(outMatrix)
	colnames(outDataFrame) = c("Hubs", "Num Interactors", "Interactors Found", "P-value", "Adjusted P-value")
	rownames(outDataFrame) = NULL
	return(outDataFrame)
	################## End Changes ###################
	
}

#######################################################################
## Save the CC for each hub-interactor pair for each condition
##
## Input
##	inHub: Name of hub
##	inInteractors: Names of interactors
##	corMatrix: Matrix of CC for all conditions
##	oFile: Output file name
##	getGeneSymb: TRUE/FALSE. If TRUE, convert Entrez IDs to gene symbols
##	isMicro: TRUE/FALSE. If TRUE, then the regulator corresponds to microRNAs
##
## Output
##	Output file
#######################################################################

############# Original code #############
# saveCorValues = function(inHub, inInteractors, corMatrix, oFile, getGeneSymb, isMicro) {
########################################
saveCorValues = function(inHub, inInteractors, corMatrix, oFile, getGeneSymb, isMicro, sampleGroups) {
########################################
	currentHub = inHub
	currentInteractors = inInteractors

	if(getGeneSymb) {
		if(!isMicro) currentHub = entrezToGeneSymbol(currentHub)
		currentInteractors = entrezToGeneSymbol(currentInteractors)
	}

	outMatrix = matrix(rep(currentHub, length(currentInteractors)), ncol=1)
	outMatrix = cbind(outMatrix, currentInteractors)
	outMatrix = cbind(outMatrix, corMatrix)
	
	if(ncol(corMatrix) == 2) {outMatrix = cbind(outMatrix, abs(corMatrix[, 1] - corMatrix[, 2]))
	
	
	
	################ Original Code ####################	
	# write(t(outMatrix), oFile, sep="\t", ncolumns=ncol(outMatrix), append=TRUE)
	##################################################
	outDataFrame = as.data.frame(outMatrix)
	rownames(outDataFrame) = NULL
	colnames(outDataFrame) = c("Hub", "Interactor", sampleGroups[1], sampleGroups[2], paste(sampleGroups[1], sampleGroups[2], sep="-"))
	  }
	return(outDataFrame)
	################## End Changes ###################
	
	
}

#######################################################################
## If during the conversion of Gene Symbols to Entrez IDs some hubs/interactors
## are missing, save them
##
## Input
##	inData: Vector of gene symbols with missing Entrez IDs
##	dataType: "Expr" - Entrez IDs missing in expression dataset
##		  "PPI_Hubs"- Entrez IDs missing in PPI hubs
##		  "PPI_Int"- Entrez IDs missing in PPI interactors
##
## Output
##	Output file
#######################################################################

generateErrorOutput = function(inData, dataType) {

	if(dataType == "Expr") write(inData, "Error_Expr.txt", ncolumns=1)	
	if(dataType == "PPI_Hubs") write(inData, "Error_PPI_Hubs.txt", ncolumns=1)	
	if(dataType == "PPI_Interactors") write(t(inData), "Error_PPI_Int.txt", ncolumns=2, sep="\t")

}

#######################################################################
## Generate a common output file that combines the results of 
## multiple expression-PPI combination analyses
##
## Input
##	fileNames: Vector of file names corresponding to individual 
##		   expression-PPI datasets results
##	outFile: File name that contains the combined information
##	metaAnalysis: Type of meta-analysis to perform
##			NULL: No meta-analysis to perform. (Default)
##			Fisher: Fisher's combined analysis
##			RankProd: Rank-prod analysis
##	rankProdItr: Number of iteration to perform for significance 
##		     of rank prod values. Default: NULL. The default 
##		     value should be used only if rank prod analysis 
##		     is not to be performed
##
## Output
##	Summarized Output file with optional columns corresponding to 
##	meta-analysis
#######################################################################

summarizeHubData = function(fileNames=NULL, outFile, metaAnalysis=NULL, rankProdItr=NULL) {

	currentNames = fileNames
	if(is.null(currentNames)) currentNames = dir(pattern="*.txt")
	
	numFiles = length(currentNames)

	allHubs = NULL
	for(i in 1:numFiles) {
	
		tempData = as.matrix(read.table(currentNames[i], header=TRUE, sep="\t"))
		allHubs = union(allHubs, tempData[, 1])
	
	} ## All files have been considered
	
	summaryMatrix = matrix(0, nrow=length(allHubs), ncol=2*length(currentNames))
	
	for(i in 1:numFiles) {
	
		indexVect = (i-1)*2
		indexVect = c(indexVect+1, indexVect+2)
	
		tempData = as.matrix(read.table(currentNames[i], header=TRUE, sep="\t"))
		allHubIndexes = match(allHubs, tempData[, 1])
		summaryMatrix[, indexVect] = tempData[allHubIndexes, c(3, 4)]		## Unadjusted p-values are saved
		
	} ## All hubs have been considered
	
	summaryMatrix = cbind(allHubs, summaryMatrix)
	
	nonExtNames = gsub(".txt", "", currentNames)
	headerInfo = c(nonExtNames[1], " ")
	for(i in 2:numFiles) headerInfo = c(headerInfo, c(nonExtNames[i], " "))
	
	headerInfo1 = c("Hub", headerInfo)
	headerInfo2 = c(" ", rep(c("Interactors Found", "P-value"), numFiles))
	
#	return(summaryMatrix)
	
	## Perform meta-analysis
	if(is.null(metaAnalysis)) {
	
		write(headerInfo1, outFile, ncolumns=length(headerInfo1), sep="\t")
		write(headerInfo2, outFile, ncolumns=length(headerInfo1), sep="\t", append=TRUE)
		write(t(summaryMatrix), outFile, ncolumns=ncol(summaryMatrix), sep="\t", append=TRUE)
	
		return(0)
	}
	
	if(metaAnalysis == "Fisher") performFisherCombinedTest(summaryMatrix, headerInfo1, headerInfo2, outFile)
	if(metaAnalysis == "RankProd") performRankProd(summaryMatrix, headerInfo1, headerInfo2, outFile, rankProdItr)
		
}

#######################################################################
## Perform Fisher's combined test
##
## Input
##	inMatrix: Matrix of p-values
##	headerInfo1: Header data to be used as first row of the o/p file
## 	headerInfo2:  Header data to be used as second row of the o/p file
## 	outFile: Output file name
##
## Output
##	Summarized output file with the last column corresponding to 
##	meta-analysis
#######################################################################

performFisherCombinedTest = function(inMatrix, headerInfo1, headerInfo2, outFile) {

	numConditions = (ncol(inMatrix) - 1)/2
	numTestIndexes = seq(from=3, to=ncol(inMatrix), by=2)
	chiProbVector = NULL

	## Perform Fisher's combined test for all hubs
	for(i in 1:nrow(inMatrix)) {
	
		exprVector = inMatrix[i, numTestIndexes]
		naCond = which(is.na(exprVector))
		currentConditions = numConditions
		
		if(length(naCond) > 0) {		
			currentConditions = currentConditions - length(naCond)
			exprVector = exprVector[-naCond]
		}
	
		exprVector = as.numeric(exprVector)
		currentProb = pchisq(q=-2*sum(log(exprVector)), df=2*currentConditions, lower.tail=FALSE)
		
		chiProbVector = c(chiProbVector, currentProb)	
	
	} ## Combined p-value calculated for all hubs
	
	## Save the Fisher's combined analysis results	
	revInfo = c(headerInfo2, "Fisher_Combined_P")	
	srcMatrix = cbind(inMatrix, chiProbVector)
	
	write(headerInfo1, outFile, ncolumns=length(headerInfo1), sep="\t")
	write(revInfo, outFile, ncolumns=length(revInfo), sep="\t", append=TRUE)
	write(t(srcMatrix), outFile, ncolumns=ncol(srcMatrix), sep="\t", append=TRUE)
	
}

#######################################################################
## Perform Fisher's combined test
##
## Input
##	inMatrix: Matrix of p-values
##	headerInfo1: Header data to be used as first row of the o/p file
## 	headerInfo2:  Header data to be used as second row of the o/p file
## 	outFile: Output file name
##
## Output
##	Summarized output file with two columns corresponding to 
##	meta-analysis
#######################################################################

performRankProd = function(inMatrix, headerInfo1, headerInfo2, outFile, rankProdItr) {

#	require(MADAM)

	## Run RankProd in parallel
	numCores = 4
	if(numCores > detectCores()) numCores = detectCores()
	cl = makeCluster(numCores)
	
	## Obtain the ranks of all hubs per experiment	
	## Set NA to the lowest rank, i.e. nrow(inMatrix)

	lowestRank = nrow(inMatrix)
	numConditions = (ncol(inMatrix) - 1)/2
	numTestIndexes = seq(from=3, to=ncol(inMatrix), by=2)
	
	rankMatrix = matrix(0, nrow=lowestRank, ncol=numConditions)
	
	for(i in 1:numConditions) {
	
		exprVector = as.numeric(inMatrix[, numTestIndexes[i]])
		rankMatrix[, i] = rank(exprVector, ties.method="min")
		
		naVector = which(is.na(exprVector))
		if(length(naVector) > 0) rankMatrix[naVector, i] = lowestRank
	
	} ## All experiments have been considered

	rankProdResults = calculateRankProduct(rankMatrix, B=rankProdItr, cluster=cl)	
	stopCluster(cl)

	## Save the RankProd results
	revInfo = c(headerInfo2, "Rank_Product_P_Value", "Rank_Product_FDR_Value")	
	srcMatrix = cbind(inMatrix, rankProdResults$p.value, rankProdResults$q.value)
	
	write(headerInfo1, outFile, ncolumns=length(headerInfo1), sep="\t")
	write(revInfo, outFile, ncolumns=length(revInfo), sep="\t", append=TRUE)
	write(t(srcMatrix), outFile, ncolumns=ncol(srcMatrix), sep="\t", append=TRUE)

}

####################################################################
## Visualize the change in correlation for the interactors of a 
## given hub
##
## Input
##	inputFile: File containing the hub-interactor pairs and 
##		   their correlation values in each condition
##	inputHub: Name of the hub
##	paletteVector: A vector of three colors
##	numThresh: Number of interactors to plot. Default is NULL
##		   which implies that all the interactors are plotted
##	condVector: Vector of the two states to plot. Default is NULL
##
## Output
##	An R graph
####################################################################

visualizeNetwork = function(inputFile, inputHub, paletteVector, numThresh=NULL, condVector=NULL) {

	inputData = read.table(inputFile, header=TRUE, sep="\t")
	
	## Determine if the number of states is greater than two
	numConditions = 2
	
	state1 = colnames(inputData)[3]
	state2 = colnames(inputData)[4]
	state3 = colnames(inputData)[5]
	
	## The minus sign is read as "."
	if(state3 != paste(state1, state2, sep=".")) numConditions = -1
	
	if(numConditions == -1) {
	
		if(is.null(condVector)) {
			print("Specify the names of the two states to be visualized")
			return(0)
		}	
		
		if(length(condVector) > 2) {
			print("Specify the names of only two states")
			return(0)
		}
		
		condIndexes = match(condVector, colnames(inputData))
		inputData = inputData[, c(1, 2, condIndexes)]
		currentColNames = colnames(inputData)
		
		corResults1 = as.numeric(inputData[, 3])
		corResults2 = as.numeric(inputData[, 4])
		inputData = cbind(inputData, corResults1 - corResults2)
		colnames(inputData) = c(currentColNames, paste(condVector[1], condVector[2], sep="-"))
	}
		
	conditionVector = colnames(inputData)[3:4]
	
	hubIndexes = grep(inputHub, inputData[, 1])
	numInteractors = length(hubIndexes)

	interactors = as.character(inputData[hubIndexes, 2])
	coeffVector1 = as.numeric(inputData[hubIndexes, 3])
	coeffVector2 = as.numeric(inputData[hubIndexes, 4])

	## If the number of interactors is more than 20, focus on the top "numThresh" interactors	
	if(numInteractors > 20 && !is.null(numThresh)) {
	
		diffInteractors = abs(as.numeric(inputData[hubIndexes, 5]))
		interactorThresh = sort(diffInteractors, decreasing=TRUE)[numThresh]
		topInteractors = which(diffInteractors >= interactorThresh)

		interactors = interactors[topInteractors]
		coeffVector1 = coeffVector1[topInteractors]
		coeffVector2 = coeffVector2[topInteractors]
				
		numInteractors = length(interactors)
	}	
		
	nwTable = matrix(rep(inputHub, numInteractors), ncol=1)
	nwTable = cbind(nwTable, interactors)
	vertexLabels = data.frame(c(interactors, inputHub))
	nwGraph = graph.data.frame(nwTable, directed=FALSE, vertices=vertexLabels)

	colorVals = getColors(paletteVector, coeffVector1, coeffVector2)
	colors1 = colorVals$actColors[1:numInteractors]
	colors2 =  colorVals$actColors[1:numInteractors + numInteractors]

	## Layout for the plots
	twoCondMatrix = matrix(1:4, nrow=2, byrow=TRUE)
	layout(twoCondMatrix, heights=c(1, 0.2))		## height (row 1) = 1, height (row 2) = 0.2
	par(mar = c(0, 2, 5, 2))
	
	## Layout for the network
	inputLayout = layout.auto(nwGraph)
	vertexColor = c(rep("lightgrey", numInteractors), "slategrey")
	
	plot(nwGraph, rescale = TRUE, layout = inputLayout
		, edge.width = 2
		, edge.color = colors1
		, vertex.frame.color = vertexColor
		, vertex.color = vertexColor
		, vertex.label.dist = 0
		, vertex.label.color = "black"
		, vertex.label = as.character(vertexLabels[, 1])
		, main = paste("Hub and interactors - ", conditionVector[1], sep=""))
	
	plot(nwGraph, rescale = TRUE, layout = inputLayout
		, edge.width = 2
		, edge.color = colors2
		, vertex.frame.color = vertexColor
		, vertex.color = vertexColor
		, vertex.label.dist = 0
		, vertex.label.color = "black"
		, vertex.label = as.character(vertexLabels[, 1])
		, main = paste("Hub and interactors - ", conditionVector[2], sep=""))	

	## Generate the spectrum chart
	numUniqueColors = length(colorVals$palColors)
	dataPoints = seq(from=1, to=numUniqueColors)

	par(mar = c(6, 3, 0, 3))
	image(matrix(dataPoints, numUniqueColors, 1), col=colorVals$palColors, axes = FALSE)
	axis(side=1, at=0, labels="Neg")
	axis(side=1, at=0.5, labels="No Cor")
	axis(side=1, at=1, labels="Pos")

}

####################################################################
## Generate the color palette using the R function maPalette
##
## Input
##	paletteVector: A vector of three colors
##	inputCoeff1: Correlation coefficients in condition 1
##	inputCoeff2: Correlation coefficients in condition 2
##
## Output
##	A list containing two elements and each element is a vector 
##	of colors
####################################################################

getColors = function(paletteVector, inputCoeff1, inputCoeff2) {

	inputCoeff = c(inputCoeff1, inputCoeff2)
	uniqueVals = sort(unique(inputCoeff))
	
	paletteColors = maPalette(low=paletteVector[1], mid=paletteVector[2], high=paletteVector[3], k=length(uniqueVals))
	actualColors = rep("white", length(inputCoeff))
	
	for(i in 1:length(uniqueVals)) {
		tempIndex = which(inputCoeff == uniqueVals[i])
		actualColors[tempIndex] = paletteColors[i]
	}
	
	return(list(actColors = actualColors, palColors = paletteColors))
				
}

####################################################################
## Filter the output hub-interactor file to consider a subset
##
## Input
##	filePrefix: File name prefix corresponding to the enriched
##		    hubs and hub-interactor pairs
##	useAdjustedProb: TRUE -> Consider adjusted p-values
##			 FALSE -> Consider unadjusted p-values
##	probThresh: p-value cut-off
##	hubNames: Vector of hub names. Default: NULL
##
## Output
##	A subset of the file corresponding to <filePrefix>_Cor.txt
####################################################################

obtainPairSubset = function(filePrefix, useAdjustedProb, probThresh, hubNames=NULL) {
	
	probColumn = 4
	if(useAdjustedProb) probColumn = 5
	
	hubResultsFile = paste(filePrefix, ".txt", sep="")
	hubPairsFile = paste(filePrefix, "_Cor.txt", sep="")
	
	resultDataHubNames = NULL
	
	if(!is.null(hubNames)) resultDataHubNames = hubNames
	else {	
		resultDataHubs = as.matrix(read.table(hubResultsFile, header=TRUE, sep="\t"))
		resultDataProb = as.numeric(resultDataHubs[, probColumn])
		resultDataHubNames = resultDataHubs[which(resultDataProb < probThresh), 1]
	}	
	
	stopifnot(length(resultDataHubNames) > 0)

	outputHeader = as.matrix(read.table(hubPairsFile, header=FALSE, sep="\t", nrows=1))
	resultDataPairs = read.table(hubPairsFile, header=TRUE, sep="\t")
	pairHubNames = as.character(resultDataPairs[, 1])
	
	outMatrix = matrix("X", nrow=1, ncol=5)
	
	## For each hub of interest, sort the correlation coefficients based on the first state
	for(hubIndex in 1:length(resultDataHubNames)) {
	
		currentInteractorIndexes = which(pairHubNames %in% resultDataHubNames[hubIndex])
		firstStateCor = as.numeric(resultDataPairs[currentInteractorIndexes, 3])	
		rankOrder = rank(firstStateCor, ties.method="random")
		orderIndexes = currentInteractorIndexes[match(1:length(rankOrder), rankOrder)]
	
		## This step is needed to ensure that the outputs in "Cor.txt" and "Cor_Signif.txt"
		## files are exactly the same in format
		tempMatrix = as.character(resultDataPairs[orderIndexes, 1])
		for(i in 2:ncol(resultDataPairs)) tempMatrix = cbind(tempMatrix, as.character(resultDataPairs[orderIndexes, i]))
		
		outMatrix = rbind(outMatrix, tempMatrix)

	}
		
	outMatrix = outMatrix[-1, ]	
	outputFileName = paste(filePrefix, "_Cor_Signif.txt", sep="")
	write(outputHeader, outputFileName, sep="\t", ncolumns=ncol(resultDataPairs))	
	write(t(outMatrix), outputFileName, sep="\t", ncolumns=ncol(resultDataPairs), append=TRUE)

}

