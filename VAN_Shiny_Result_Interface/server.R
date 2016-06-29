
## Extra functions for better network visualisations in Shiny, using iGraph. 
source("VAN_Plotting.R")




shinyServer(function(input, output, session) {

        readingExistingFiles  = eventReactive(
          eventExpr = input$existingButton, 
          valueExpr = {
            if(is.null(input$exsitingProbFile) || is.null(input$exsitingCorFile)) {return(NULL)}
            else{probFile = read.delim(input$exsitingProbFile$datapath, header=T)
            corFile = read.delim(input$exsitingCorFile$datapath, header=T)
            twoFiles = list(probFile = probFile,corFile = corFile)
            return(twoFiles)
          }
        })
    
    output$probResult = renderDataTable({
      return(readingExistingFiles()$probFile)
      })
    
    output$corResult = renderDataTable({
      return(readingExistingFiles()$corFile)
      })
    
######### Plotting functionalities ###############
####################################################################################
    
    observe({
      QualifiedHubs = as.character((readingExistingFiles()$probFile)[,1])  
      updateSelectInput(session, "InputPlotHub", choices = QualifiedHubs)
      updateSliderInput(session, "nTopGenes", max = length(QualifiedHubs))
    })
    
    
    
    
    # This plots the first condition in the data.   
    output$GraphCond1 = renderPlot({
      corResult = readingExistingFiles()$corFile
      
      hubNetwork = subset(corResult, Hub == input$InputPlotHub)
      
      # These two correlation measures has a sign attached. Hence the use of colours in plotting the networks
      if (input$CorType == "PCC" || input$CorType == "Spearman"){
        return(plotBWR(cornetwork = hubNetwork, condition = colnames(hubNetwork)[3]))
      }
      
      # These two correlation measures does not have a sign attached. Hence the use of Black/White in plotting the networks
      if (input$CorType == "Distance" || input$CorType == "MIC"){
        return(plotBW(cornetwork = hubNetwork, condition = colnames(hubNetwork)[3]))
      } 
    })
    
    
    # This plots the second condition in the data.  
    output$GraphCond2 = renderPlot({
      corResult = readingExistingFiles()$corFile
      
      
      hubNetwork = subset(corResult, Hub == input$InputPlotHub)
      
      if (input$CorType == "PCC" || input$CorType == "Spearman"){
        return(plotBWR(cornetwork = hubNetwork, condition = colnames(hubNetwork)[4]))
      }
      
      
      if (input$CorType == "Distance" || input$CorType == "MIC"){
        return(plotBW(cornetwork = hubNetwork, condition = colnames(hubNetwork)[4]))
      } 
    })
    
    
    # This plots the abs(first condition - second condition) in the data.
    # Since the difference is always non-negative, we use the black-white colouring again. 
    output$GraphCondDiff = renderPlot({
      corResult = readingExistingFiles()$corFile
      
      hubNetwork = subset(corResult, Hub == input$InputPlotHub)
      return(plotBW.Diff(cornetwork = hubNetwork, condition = colnames(hubNetwork)[5]))
    })
    
    
    output$rcytoscapejsPlot  = renderRcytoscapejs({
      corResult = readingExistingFiles()$corFile
      if (nrow(corResult)>2000){
        return(NULL) 
      } else {
        network = plot.rcytoscapejs(corResult = corResult)
        return(rcytoscapejs(network$nodes, network$edges, showPanzoom=TRUE))    
        }
    })
    
    
    
    output$overallGraph = renderPlot({
      corResult = readingExistingFiles()$corFile
      set.seed(1)
      # This forms an unweighted, undirected network
      graph = graph.data.frame(corResult[,c("Hub","Interactor")], directed=F) 
      
      result = plot.igraph(graph, vertex.size = 1, 
                  vertex.label = NA,
                  vertex.color = "#8f98ff")
      return(result)
    })
    
    
    
    
    
    output$scatterPlot = renderScatterD3({
      probResult = readingExistingFiles()$probFile
      plotProbResult = probResult[,c("Hubs","Interactors.Found","P.value")]
      plotProbResult = plotProbResult[order(plotProbResult[,"P.value"]),]
      
      nTopGenes = input$nTopGenes
      
      return(scatterD3(
        x = plotProbResult[1:nTopGenes,"Interactors.Found"],
        y = plotProbResult[1:nTopGenes,"P.value"],
        lab = plotProbResult[1:nTopGenes,"Hubs"],
        col_var = ifelse(plotProbResult[1:nTopGenes,"P.value"]<=0.05,"p<=0.05","p>0.05"),
        colors = c("#FF0000","#0080FF"),
        point_opacity = 0.6,
        labels_size = 12
      ))
    }) # End the scatterPlot
  
    
    output$HubNetwork = renderPlot({
      corResult = readingExistingFiles()$corFile
      probResult = readingExistingFiles()$probFile
      return(plot.HubNetwork(corResult = corResult, probResult = probResult))
      })
    
    
    
}) # End shinyServer
