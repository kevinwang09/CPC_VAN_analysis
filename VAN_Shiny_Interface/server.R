options(shiny.maxRequestSize = 20*1024^2)


## The four main R scripts in VAN
source("InteractomeMultiThread_Func.R")
source("ReadWrite_Func.R")
source("RankMethods.R")
source("Metadata_Func.R")

## Extra functions for better network visualisations in Shiny, using iGraph. 
source("VAN_Plotting.R")




shinyServer(function(input, output, session) {
  ## This function checks if all the necessary inputs are present. 
  ## Directly goes to the main function computeVAN.
  
  
  
    CheckAllInputs = function(){
      # If any of the inputs are empty, we will halt the computation at some point. 
      ALLCLEAR = TRUE
      if (length(input$InputHub)==0) ALLCLEAR = FALSE
      if (is.null(input$ExpressionFile)==TRUE) ALLCLEAR = FALSE
      if (is.null(input$PPIFile)==TRUE) ALLCLEAR = FALSE
      if (length(input$HubSize)==0) ALLCLEAR = FALSE
      return(ALLCLEAR)
    }
    
    
    
    # Observing the input hubs
    observe({
      if(is.null(input$PPIFile)==TRUE) {
        QualifiedHubs = "Please upload PPI first." 
        # Prior to uploading PPI map, a message to the user 
        updateSelectInput(session, "InputHub", choices = QualifiedHubs)
      } else {
        CompleteHubs = as.character(read.table(input$PPIFile$datapath, header=T)[,1])
        QualifiedHubs = names(which(table(CompleteHubs) >= input$HubSize)) # The HubSize is used to do a very quick filtering job.
        # Once the PPImap is uploaded, we quickly find out the name of those hubs, and ask
        # the user what hubs do they want to compute the significance for. 
        QualifiedHubs = c("All hubs", sort(unique(QualifiedHubs)))
        # Computing all hubs is also a valid option.
        updateSelectInput(session, "InputHub", choices = QualifiedHubs)
      }
    })
    
    
#     
#     existing = reactiveValues()
#     existing$value = FALSE
# 
#     observeEvent(eventExpr = input$existingButton, handlerExpr = {
#       existing$value  = TRUE
#     })
#     
#     observeEvent(eventExpr = input$expressionButton, handlerExpr = {
#       existing$value  = FALSE
#     })
#     
#     

    
    
    ComputingVAN = eventReactive(eventExpr = input$expressionButton, valueExpr = {
      # This part is only activated via clicking the button. 
      if (CheckAllInputs() == FALSE){ stop("Enter all fields to start") } 
      # If any inputs are missing, then we display a warning and do not proceed further. 
      else{
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...', value = 0, {
                       ## Buried deep inside the InteractomeMultiThread_Func.R
                       # there is a line of incProgress(...) that works with this function
                       # The idea is to tell the user the progress of execution
                       # in terms of how many hubs are left.
                       
                       if ("All hubs" %in% input$InputHub) {hubVect = NULL} # If All hubs were selected, we hubVect of VAN's input as NULL
                       else {hubVect = as.character(input$InputHub)}
                       
                       Results = identifySignificantHubs(exprFile = input$ExpressionFile$datapath,
                                                         mapFile = input$PPIFile$datapath,
                                                         outFile = "Testing_Short_PPI.txt", 
                                                         labelIndex = 1, 
                                                         randomizeCount = input$NumPermTests,
                                                         hubSize = input$HubSize,
                                                         assocType = input$CorType,
                                                         hubVect = hubVect,
                                                         RunShinyAPP = TRUE)
                       # Standard VAN inputs. Except the RunShinyAPP variable. 
                     })
        closeAllConnections()
        return(Results)
      }
      
    })# End eventReactive and valueExpr
##################################################################################################    

    
    
    # This probResult lists out all the hubs and their corresponding P-value via permutation test. 
    output$probResult = renderDataTable({
      return(ComputingVAN()$probResult)
    })
    
    
    
    # This corResult lists out pairwise correlations. 
    output$corResult = renderDataTable({
      return(ComputingVAN()$corResult)
    })
    
    
    
    
  
######### Plotting functionalities ###############
####################################################################################
    
    observe({
      QualifiedHubs = as.character((ComputingVAN()$corResult)[,1])  
      updateSelectInput(session, "InputPlotHub", choices = QualifiedHubs)
    })
    
    
    
    
    # This plots the first condition in the data.   
    output$GraphCond1 = renderPlot({
      corResult = ComputingVAN()$corResult
      
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
      corResult = ComputingVAN()$corResult

      
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
        corResult = ComputingVAN()$corResult
      
      hubNetwork = subset(corResult, Hub == input$InputPlotHub)
      return(plotBW.Diff(cornetwork = hubNetwork, condition = colnames(hubNetwork)[5]))
    })
    
    
    
    
    output$rcytoscapejsPlot  = renderRcytoscapejs({
      corResult = ComputingVAN()$corResult
      network = plot.rcytoscapejs(corResult = corResult)
      return(rcytoscapejs(network$nodes, network$edges, showPanzoom=TRUE))
      
    })
    
    
  
  
  
}) # End shinyServer
