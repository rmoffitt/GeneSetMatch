#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#shiny::runApp(appDir = "~/ODIS2/Shiny", display.mode="showcase")
library(shiny)

source("~/ODIS2/Shiny/prepareLLB.R")
LLB <- prepareLLB() 
full_Snyder <- prepareSnyder()

options(shiny.maxRequestSize = 30*1024^2) #30MB
completed = FALSE
analysisres <- NULL
GSEAres <- NULL
dataset <- LLB





# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("ODIS"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          # Selecting a data set
            h3("Upload dataset or use example dataset"),
            fileInput("dataFile", #name of item in app
                      "File input", #name appearing near button
                      multiple = FALSE,
                      accept = c(".Rdata", ".csv", "text/csv", "plain"),
                      placeholder = "select expression file"),
            
            radioButtons("dataset", "select dataset",
                         choices = list("Liver/Lung/Brain" = "LLB", "Snyder" = "snyder", "Null" = "3"),
                         selected = "LLB"),
            
            #--------------------------------------------------------------
            # selecting an analysis
            h3("Select Analysis"),
            radioButtons("analysisChoice", "Analysis",
                         choices = list("DESeq" = "deseq", "NMF" = "nmf", "FindAllMarkers" = "fam",
                                        "etc" = "null"),
                         selected = "deseq"),
            # Create options based on method selection
            conditionalPanel("input.analysisChoice == 'deseq'",
                             p("deseq options go here"),
                             sliderInput("placeholder", "place holder:",
                                         min = 2, max = 20, #change max to ncol? Just cap at something reasonable?
                                         value = 3,
                                         ticks = F
                             )
            ),
            conditionalPanel("input.analysisChoice == 'nmf'",
                             p("You can test several ranks, resulting in nmf rank plots, or simply
                                select a rank if you already know what to use."),
                             sliderInput("rankToTest", "Ranks to test:",
                                         min = 2, max = 20, #change max to ncol? Just cap at something reasonable?
                                         value = c(3, 10),
                                         ticks = F
                             ),
                             actionButton("rankTest", "test different ranks"),
                             br(),
                             br(),
                             p("-----------------------"),
                             sliderInput("nmfRank", "rank to use:", #issue -- ideally would not display blue on color bar, but fine for now
                                         min = 2, max = 20,
                                         value = 3,
                                         ticks = F
                             )
            ),
            conditionalPanel("input.analysisChoice == 'fam'",
                             p("FindAllMarkers options go here"),
                             sliderInput("placeholder", "place holder:",
                                         min = 2, max = 20, #change max to ncol? Just cap at something reasonable?
                                         value = 3,
                                         ticks = F
                             )
            ),
            actionButton("runAnalysis", "run analysis")
        ),

        # Tabsets with output from analysis, ideally will only appear after the data has been generated?
        mainPanel(
          tabsetPanel(id = "mainTabs",
                      type = "tabs",
                      tabPanel("Preview raw data", DT::dataTableOutput("rawData")),
                      tabPanel("analysis output", plotOutput(("analysisPlot")), uiOutput("GSEAbutton")), #updateTabsetPanel to change name?
                      tabPanel("GSEA output", uiOutput("genesetRadios"), DT::dataTableOutput("gseaTable"), actionButton("makeHeatmap", "generate heat map")),
                      tabPanel("ODIS Heatmap")
          )
        )
    )
)



# ============================================================================================================
# ============================================================================================================


# Define server logic 
server <- function(input, output, session) {
    
  # https://stackoverflow.com/questions/52385295/get-the-name-of-an-uploaded-file-in-shiny
  #output$fileName <- renderText({
  #    # Test if file is selected
  #    if (!is.null(input$dataFile$datapath)) {
  #        # Extract file name (additionally remove file extension using sub)
  #        return(basename(input$dataFile$name))
  #    } else {
  #        return(NULL)
  #    }
  #})
  
  output$analysisName <- renderText({
    return(input$analysisChoice)
  })

  
  # dataset selection with input dataset
  observe({
    if (is.null(input$dataFile))
      return()
    e = new.env()
    name <- load(input$dataFile$datapath, envir = e)
    inputfile <<- e[[name]]

    # Can also set the label and select items
    updateRadioButtons(session, "dataset",
                       label = paste("Upload dataset or use example dataset"),
                       #choices = list("Liver Lung Brain" = 1, "Choice 2" = 2, input$dataFile$name = 3),
                       choiceNames = c("Liver/Lung/Brain", "Snyder", "null", input$dataFile$name),
                       choiceValues = c("LLB", "snyder", "3", "file"),
                       selected = "file")
  })
  
  # set which dataset is being used as a global variable
  observe({ # issue???????!???
    dataset <<- switch(input$dataset,
                       "LLB" = LLB,
                       "snyder" = full_Snyder,
                       "3" = NULL, # issue -- delete or add data set
                       #"4" = read.csv(input$dataFile$datapath)),
                       "file" = inputfile)
    print(input$dataset)
  })

  # issue -- Would like this to save the desired data table as a variable to be accessed later
  output$rawData <- DT::renderDataTable(
    #dataset, # not reactive.
    switch(input$dataset,
           "LLB" = LLB,
           "snyder" = full_Snyder$ex,
           "3" = NULL, # issue -- delete or add data set
           #"4" = read.csv(input$dataFile$datapath)),
           "file" = inputfile),
    options = list(scrollX = TRUE)
  )
    
  # ---------------------------------------------------------------
  # ---------------------- ANALYSIS -------------------------------
  # ---------------------------------------------------------------
    
  #NMF rank test
  observeEvent(input$rankTest, {
    #dataset <- head(dataset, 100)
      
      
    showNotification("Running nmfEstimateRank. This will take a minute.")
    ranks <- seq(input$rankToTest[1], input$rankToTest[2])
    #ranks <- input$rankToTest[1]:input$rankToTest[2]
    print(ranks)
    nmfTest <- estimateRankShiny(dataset, ranks)
    #nmfTest <- NMFforShiny(dataset, ranks, nrun = 10)
    print("done")
    output$analysisPlot <- renderPlot(plot(nmfTest))
    updateTabsetPanel(session, "mainTabs", selected = "analysis output")
      
  })
    
  # run analysis based on selection
  observeEvent(input$runAnalysis, {
    showNotification("Running analysis. This will take a minute.")
    #print(dataset)
    #dataset <- LLB    # issue ?
    completed <- FALSE
      
    if(input$analysisChoice == "nmf"){
      print("running NMF pipeline")
      source("~/ODIS2/Shiny/shinyNMF.R") # issue
      analysisres <<- NMFforShiny(dataset, input$nmfRank)

      # issue -- may no longer be needed by use of sliders
      if(class(analysisres) == "character"){
        showNotification(analysisres)
      }
      else{
        # issue -- determine how to format res
        output$analysisPlot <- renderPlot(consensusmap(analysisres))
        updateTabsetPanel(session, "mainTabs", selected = "analysis output")
        completed <- TRUE
      }
      print("NMF completed")
      #View(res)
    } 
    # End of NMF analysis
    else{
      showNotification("not yet implemented")
      # issue -- TBD
    }
      
    if(completed){
      output$GSEAbutton <- renderUI({
        tagList(checkboxGroupInput("pathwaySelection", "Select pathways for GSEA",
                           choices = c("C1", "C2", "C3","C4","C5","C6","C7","C8", "custom (this is a placeholder and selecting it breaks the code)"),
                           selected = NULL),
        actionButton("runGSEA", "run GSEA"))
      })
    }
  })
    
  
  # -----------------------------
  # --------- GSEA --------------
  # -----------------------------
  observeEvent(input$runGSEA, {
    showNotification("running GSEA... functioning!? :0") # :(
    print(input$pathwaySelection)
    
    source("~/ODIS2/Shiny/GSEAfunctions.R") # issue
    species <- switch (input$dataset,
      "LLB" = "Rattus norvegicus",
      "snyder" = "Mus musculus",
      "3" = NULL,
      "file" = "Mus musculus" # issue -- need to create some sort of input widget
    )
    print(species)
    genesets <- prepareGenesets(species, input$pathwaySelection) # issue -- allow user to select genesets. For now just testing on C8...
    #View(res)
    gene_list <- sort(analysisres@fit@W[,1], decreasing = TRUE) # issue -- place holder to test with NMF, need to determine res structure
    #print(gene_list)
    GSEAres <- ODISGSEA_helper(gene_list = gene_list, gmtList = genesets, pval = 1) #maintains listed structure for compatibility with ODISHEATMAP, but should only ever have one list item
    #View(GSEAres)
    
    output$genesetRadios <- renderUI({
      radioButtons("viewGeneset", "Geneset to view",
                   choiceNames = names(GSEAres[[1]]),
                   choiceValues = names(GSEAres[[1]]),
                   select = names(GSEAres[[1]])[1]
                   )
    })
    print("generating table")
    # issue "Warning: Error in [[: attempt to select less than one element in
    # get1index" I believe this may due to input$viewGeneset not existing then
    # datatable is generated, then it does not continue to show the error?
    output$gseaTable <- DT::renderDataTable(GSEAres[[1]][[input$viewGeneset]][["Results"]],
                                            options = list(scrollX = TRUE)) # issue -- make this viewable, and then nicely navigable 
    updateTabsetPanel(session, "mainTabs", selected = "GSEA output")
  })
  
  # ------------------------------------
  # ----------- Heatmap ----------------
  # ------------------------------------
  observeEvent(input$makeHeatmap, {
    showNotification("You clicked the button! :)")
  })
  
    
}

# Run the application 
shinyApp(ui = ui, server = server)
