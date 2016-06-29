shinyUI(
  fluidPage(
  titlePanel("VAN --- Differential Variability Analysis via Network Approach"),
  h3("A Shiny Interface for computation on expression files"),
  
  sidebarLayout( # The first component has a sidebarLayout
  sidebarPanel(
    h3("File inputs"),
    hr(),
    actionButton("expressionButton", label = "Click here to start VAN"),
    hr(),
    
    fileInput("ExpressionFile", label = h4("Expression file input")),
    br(),
    
    fileInput("PPIFile", label = h4("PPI file input")),
    "-----------------------------------------",
    
    h3("VAN computation parameters"),
    
    numericInput(inputId = "HubSize", label = h4("Minimum number of interactors to qualify as a hub"), value = 5, min = 1),
    
    selectInput("InputHub", label = h4("Select hub gene(s) of interest"), 
                choices = "QualifiedHubs", 
                multiple = TRUE),
    
    radioButtons("CorType", label = h4("Type of correlation"),
                 choices = list(
                   # "Taylor Correlation" = "TCC", 
                   "Pearson Correlation" = "PCC",
                   "Spearman Correlation" = "Spearman",
                   "Distance Correlation" = "Distance",
                   "Maximal Information Criterion" = "MIC"),
                 selected = "PCC"),
    
    numericInput(inputId = "NumPermTests", label = h4("Number of permutation test iterations"), value = 1000, min = 10),
    
    "-----------------------------------------",
    
    h3("VAN plotting parameters"),
    
    selectInput("InputPlotHub", label = h4("Select only one hub to plot"), 
                choices = "QualifiedPlotHubs", 
                multiple = FALSE)
  ), # End sidebarPanel
####### For the first component, moving from sidebarPanel into the mainPanel
    mainPanel(
      tabsetPanel(
        tabPanel("Summary Tables",
                 h3("Summary of hub results"),
                 dataTableOutput("probResult"),
                 
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
                 
                 h3("Interactive Grand Network"),
                 rcytoscapejsOutput("rcytoscapejsPlot", height="600px")
        ), # End second tabPanel
        tabPanel("Documentations",
                 p("Given appropriate inputs, this app acts as an interactive interface to the 'identifySignificantHubs' function in the VAN package."),
                 p("We refer users to the original VAN user guide for relevant instructions."),
                 h2("Inputs"),
                 h3("Input files"),
                  p("This app require two input files."),
                    h4("Expression file input"),
                      p("The first is the gene expression file in txt format. VAN requires a very specific formatting for this file. 
                        We refer users to the original VAN user guide for more information"),
                    h4("PPI file input"),
                      p("The second is the protein-protein interaction file, which contains a set of hub-interactor pairs To be investigated."),
                 h3("Input controls"),
                    h4("Minimum number of interactors to qualify as a hub"),
                      p("In order for a gene to be considered as a 'hub', we need to specify how many interactors it must have relative to the PPI file."),
                    h4("Select hub gene(s) of interest"),
                      p("The default of this option is to consider all hubs in the uploaded PPI file."),
                      p("However, if only some hubs are of interests, then after the PPI file has been uploaded, user can select a list of hubs manually."),
                      p("We do NOT recommended overwhelmingly large number of hubs in calculations, since this often is not needed and takes up precious computational resources."),
                    h4("Type of correlation"),  
                      p("We provided four different types of correlation measures. Users are free to choose which ever one that is suitable during computation. 
                        This also has an effect on color scheme during visualisations. See below."),
                    h4("Number of permutation test iterations"),
                      p("In order to establish the significance of a hub, we will need to perform the non-parametric permutation test on the patient conditions."),
                      p("We recommend 1000 re-samplings here. But it should be noted increasing the number of re-samplings will greatly increase the computation burden; often with limited increase in precision."),
                    h4("Select only one hub to plot"),
                      p("Users are free to select which hub to be plotted, once the computations has terminated."),
                 
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
                 )
      ) # End tabsetPanel
    ) # mainPanel
                      
  ) # End first component's sidebarLayout 

) # End fluidPage
) # End shinyUI