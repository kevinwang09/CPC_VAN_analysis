################ Plotting coloured edges, for Pearson and Spearman ###################
plotBWR = function(cornetwork, condition){
  ncol = 7
  rbPal <- colorRampPalette(c('blue','white','red'))
  fullcol = rbPal(ncol)
  # The three colours are extrapolated into 100 full spectrum of colours. 
  # fullcol is a character vector, with blue as its first element and red its last. 
  
  
  set.seed(1)
  # This forms an unweighted, undirected network
  graph = graph.data.frame(cornetwork[,c("Hub","Interactor")], directed=F) 
  
  # The weights are the numerical values correspond to some condition
  CorrVal = as.numeric(as.character(cornetwork[,as.character(condition)]))
  
  
  # The formula had been worked out to map the numerical values correspond to the conditions into the colors
  # CorrVal range: [-1,1]
  edgecolor = fullcol[round(0.5*(ncol-1)*CorrVal+0.5*(ncol+1))]
  
  
  finalplot = plot.igraph(graph, main = as.character(condition), 
                          edge.width = abs(15*CorrVal), # this was chosen for aesthetic reasons
                          vertex.size = 25, # this was chosen for aesthetic reasons
                          edge.color = edgecolor, 
                          vertex.color = "#FFF68F") # this was chosen for aesthetic reasons
  return(finalplot)
}
################################################################################







################ Plotting BW edges, For Differences and for MIC and Dcor (2007 version and absolute value of the 2013 version) ###################
plotBW = function(cornetwork, condition, CorrDiff = FALSE){
  ncol = 7
  rbPal <- colorRampPalette(c('white','black'))
  fullcol = rbPal(ncol)
  # The three colours are extrapolated into 100 full spectrum of colours. 
  # fullcol is a character vector, with blue as its first element and red its last. 
  
  
  set.seed(1)
  # This forms an unweighted, undirected network
  graph = graph.data.frame(cornetwork[,c("Hub","Interactor")], directed=F) 
  
  # The weights are the numerical values correspond to some condition
  CorrVal = as.numeric(as.character(cornetwork[,as.character(condition)]))
  
  
  # The formula had been worked out to map the numerical values correspond to the conditions into the colors
  # CorrVal range: [0,1]
  edgecolor = fullcol[round( (ncol-1)*CorrVal + 1)]
  
  # If we are plotting differences of two correlation measures, with values in [0,1]
  # Then, all edges are shown with the same width. 
  # If it is just a correlation measure on a condition, then we will weight the edges accordingly
    
  
  finalplot = plot.igraph(graph, main = as.character(condition), 
                            edge.width = 15*abs(CorrVal), # this was chosen for aesthetic reasons
                            vertex.size = 25, # this was chosen for aesthetic reasons
                            edge.color = edgecolor, 
                            vertex.color = "#FFF68F") # this was chosen for aesthetic reasons
  
  return(finalplot)
}
####################



################ Plotting BW edges, For Differences and for MIC and Dcor (2007 version and absolute value of the 2013 version) ###################
plotBW.Diff = function(cornetwork, condition){
  ncol = 7
  rbPal <- colorRampPalette(c('white','black'))
  fullcol = rbPal(ncol)
  # The three colours are extrapolated into 100 full spectrum of colours. 
  # fullcol is a character vector, with blue as its first element and red its last. 
  
  
  set.seed(1)
  # This forms an unweighted, undirected network
  graph = graph.data.frame(cornetwork[,c("Hub","Interactor")], directed=F) 
  
  # The weights are the numerical values correspond to some condition
  diffVal = as.numeric(as.character(cornetwork[,as.character(condition)]))
  
  
  # The formula had been worked out to map the numerical values correspond to the conditions into the colors
  # diffVal range: [0,m]
  # color value: [1,n]
  m = max(diffVal)
  
  edgecolor = fullcol[round(  diffVal*(ncol-m) + m  )]
  
  # If we are plotting differences of two correlation measures, with values in [0,1]
  # Then, all edges are shown with the same width. 
  # If it is just a correlation measure on a condition, then we will weight the edges accordingly
  
    finalplot = plot.igraph(graph, main = as.character(condition), 
                            edge.width = 13, # this was chosen for aesthetic reasons
                            vertex.size = 25, # this was chosen for aesthetic reasons
                            edge.color = edgecolor, 
                            vertex.color = "#FFF68F") # this was chosen for aesthetic reasons
  
  return(finalplot)
}
####################
plot.rcytoscapejs = function(corResult){
  
  allNodes = unique(c(as.character(corResult[,1]),as.character(corResult[,2])))
  nodeData = data.frame(id = allNodes, name = allNodes, stringsAsFactors = FALSE)
  sourceNodes = as.character(corResult[,1])
  targetNodes = as.character(corResult[,2])
  edgeData = data.frame(source = sourceNodes,
                        target = targetNodes,
                        stringsAsFactors = FALSE)
  
  nodeData$color <- rep("#00FF00", nrow(nodeData)) 
  nodeData$color[which(nodeData$id %in% unique(sourceNodes))] <- "#FF0000"
  
  network = createCytoscapeJsNetwork(nodeData, edgeData)
  return(network)
}
