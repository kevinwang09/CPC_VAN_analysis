shinyUI(
  fluidPage(
    titlePanel("VAN --- Differential Variability Analysis via Network Approach"),
    h3("A Shiny Interface for computation on expression files"),
    sidebarLayout( # The first component has a sidebarLayout
      sidebarPanel(
        actionButton("existingButton", label = "Click here to visualise VAN results"), 
        
        fileInput("exsitingProbFile", label = h3("Exisiting Prob File")),
        br(),
        
        fileInput("exsitingCorFile", label = h3("Exisiting Cor File")),
        hr(),
        
        "-----------------------------------------",
        
        h3("VAN computation parameters"),
      
        radioButtons("CorType", label = h4("Type of correlation"),
                     choices = list(
                       # "Taylor Correlation" = "TCC", 
                       "Pearson Correlation" = "PCC",
                       "Spearman Correlation" = "Spearman",
                       "Distance Correlation" = "Distance",
                       "Maximal Information Criterion" = "MIC"),
                     selected = "PCC"),
        

        "--------------------------------------------",
        h3("VAN plotting parameters"),
        selectInput("InputPlotHub", label = h4("Select only one hub to plot"), 
                    choices = "QualifiedPlotHubs", 
                    multiple = FALSE),
        sliderInput("nTopGenes", label = h4("Top n genes on a scatter plot (number of Interactors vs P-value)"), 
                    min = 1, max = 100, value = 10)
        
      ),
###################################################################      
      mainPanel(
        tabsetPanel(
          tabPanel("Summary Tables",
                   h3("Summary of hub results"),
                   dataTableOutput("probResult"), # Use selection of table rows. https://rstudio.github.io/DT/shiny.html
                   
                   h3("Summary of correlation results"),
                   dataTableOutput("corResult")
          ), # End first tabPanel
          ################ Plotting Functionalities ###############################    
          tabPanel("Network Plots",
                   h3("First condition network plot"),
                   plotOutput("GraphCond1",height="600px"),
                   
                   h3("Second condition network plot"),
                   plotOutput("GraphCond2",height="600px"),
                   
                   h3("Edge difference plot"),
                   plotOutput("GraphCondDiff",height="600px"),
                   
                   h3("Num. Inetactor vs P-value scatterplot"),
                   scatterD3Output("scatterPlot",height="600px"),
                   
                   h3("Overall Hubs-Network"),
                   plotOutput("HubNetwork",height="600px"),
                   
                   h3("Interactive Grand Network"),
                   rcytoscapejsOutput("rcytoscapejsPlot", height="600px")

          ), # End second tabPanel
          tabPanel("Documentations",
                   p("Given appropriate inputs, this app will produce a set of outputs from the VAN package to visualise network and hub significance."),
                   p("Users must run the 'identifySignificantHubs' function in the VAN package and upload the output files."),
                   h2("Inputs"),
                    h3("Input files"),
                      p("This app require two input files."),
                      
                      h4("Exisiting Prob File"),
                        p("The first is the probability file, which records the hubs and their significance."),
                   
                      h4("Exisiting Cor File"),
                        p("The second is the correlation file, which records all hub-interactor pairs under patient conditions and the corresponding correlation values."),
                    h3("Input controls"),

                      h4("Type of correlation"),
                        p("Depending on which correlation were used during the calculation, users should select the corresponding correlation type when plotting. 
                          The 'Type of correlation' controls the color scheme used when visualising individual hub."),

                      h4("Select only one hub to plot"),
                        p("Users are free to select which hub to be plotted, once the computations has terminated."),
                   
                      h4("Top n genes on a scatter plot (number of Interactors vs P-value)"),
                        p("Since the high significance are often biased towards the hubs with fewer interactors, 
                            we recommend users to consider using the slider bar to explore some general structure between significance and number of interactors."),
                  
                   h2("Outputs"),
                    h3("Visualisations"),
                      p("This app only support visualisations for gene expression profiles under two-conditions."),
                      
                      h4("First condition network plot and Second condition network plot"),                   
                        p("The first two plots are hub-interactor network for an user-selected hub, under two different conditions."),
                        p("The thickness and colour intensity of a particular edge are linearly proportional to the magnitude of the correlation value between the hub and the interactor."),                   
                      
                      h4("Edge difference plot"),
                        p("This plot is based on thre previous two. It is simply the absolute difference in the correlation values, with larger difference gainly a more intense color."),
                        p("For Pearson and Spearman correlations, the edges are coloured using blue-white-red color scheme to highlight the signs of the correlation values; 
                          with blue representing a negative value and red representing a positive value."),
                      
                      h4("Num. Inetactor vs P-value scatterplot"),
                        p("This plot can be used to visualise the overall hub significance and its relationship with the number of interators a hub has. 
                          This plot is fully interactive with the input 'op n genes on a scatter plot (number of Interactors vs P-value)'"),
                        p("Since the high significance are often biased towards the hubs with fewer interactors, 
                            we recommend users to consider using the slider bar to explore some general structure between significance and number of interactors."),
                      
                      h4("Overall Hubs-Network"),
                        p("This plots how every hubs are linked to each other. Size of the radius is proportional to the -log(p-value) of that hub. The hubs with p-value less than 0.05 are highlighted in red. 
                          This plot is particularly useful if the networks are small. "),
                      
                      h4("Interactive Grand Network"),
                        p("If the inputed correlation file has less than 2000 rows, i.e. less than 2000 hub-interactor pairs; then, 'Interactive Grand Network' is an interactive tool to investigate the grand network. 
                          If the number of hub-interactor pairs exceeds this, the output will be suppressed due to computational difficulties.
                          If such highly complex network visualisation is needed, we recommend the free software Cytoscape.")
                   ) # End the tabPanel on documentation.
        ) # End tabsetPanel
      ) # mainPanel
      
    )
  ) # End fluidPage
) # End shinyUI