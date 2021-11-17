#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
options(shiny.maxRequestSize = 30*1024^2) #30MB
source("~/ODIS2/inst/shiny/shinyNMF.R") # :(

# -------------------------------------------------------------
# ---------------------------- UI -----------------------------
# -------------------------------------------------------------
ui <- fluidPage(
    # Application title
    titlePanel("ODIS"),
    
    # --------------
    # -- SIDE BAR --
    # --------------
    column(3,
        wellPanel(
            # dataset UI
            h3("Upload dataset or use example dataset"),
            fileInput("dataFile", #name of item in app
                      "File input", #name appearing near button
                      multiple = FALSE,
                      accept = c(".Rdata", ".csv", ".Rds"),
                      placeholder = "select expression file"),
            
            radioButtons("datasetRadio", "select dataset",
                         choices = list("Liver/Lung/Brain" = "LLB", "Snyder" = "snyder"),
                         selected = "LLB")
        ),
        #=====================================================================================
        # analysis UI
        wellPanel(
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
        #============================================================================================
        wellPanel(
            h3("GSEA")
        )
    ),

            
    # ----------------
    # -- MAIN PANEL --
    # ----------------
    column(9,
           tabsetPanel(id = "mainTabs",
                       type = "tabs",
                       tabPanel("Preview raw data", DT::dataTableOutput("rawData")),
                       # issue -- set plot output size to prevent weird stretching
                       tabPanel("analysis output", plotOutput(("analysisPlot")), uiOutput("GSEAbutton")), #updateTabsetPanel to change name?
                       tabPanel("GSEA output", uiOutput("genesetRadios"), DT::dataTableOutput("gseaTable"), actionButton("makeHeatmap", "generate heat map")),
                       tabPanel("ODIS Heatmap")
           )
        )
    #)
)

# ------------------------------------------------------------
# ------------------------ SERVER ----------------------------
# ------------------------------------------------------------
server <- function(input, output, session) {
    # -------------
    # -- GLOBALS --
    # -------------
    globals <- reactiveValues(dataset_list = list(LLB = readRDS("~/ODIS2/data/LLB.Rds"), # issue -- build and set loads to ODIS2::dataname
                                              snyder = readRDS("~/ODIS2/data/full_snyder_ex.Rds"),
                                              file = NULL), 
                              analysisres = NULL,
                              GSEAres = NULL)
    
    # -------------
    # -- DATASET --
    # -------------
    observeEvent(input$datasetRadio, {
        print(input$datasetRadio)
    })
    
    uploadedFileName <- eventReactive(input$dataFile, {
        input$dataFile$name
    })
    
    observeEvent(input$dataFile, {
        updateRadioButtons(session, "datasetRadio",
                           label = paste("Upload dataset or use example dataset"),
                           choiceNames = c("Liver/Lung/Brain", "Snyder", uploadedFileName()),
                           choiceValues = c("LLB", "snyder", "file"),
                           selected = "file")
        
        #load in file
        # issue -- check what type of file, want to accept Rdata/Rds dataframes and csv
        e = new.env()
        globals$dataset_list$file <- e[[load(input$dataFile$datapath, envir = e)]]
        #if(csv)
        #readcsv
        })

    
    output$rawData <- DT::renderDataTable(
        globals$dataset_list[[input$datasetRadio]],
        options = list(scrollX = TRUE)
    )
    
    # --------------
    # -- ANALYSIS --
    # --------------
    
    #NMF rank test
    observeEvent(input$rankTest, {
        showNotification("Running nmfEstimateRank. This will take a minute.")
        ranks <- seq(input$rankToTest[1], input$rankToTest[2])
        print(ranks)
        nmfTest <- estimateRankShiny(globals$dataset_list[[input$datasetRadio]], ranks)
        print("done")
        output$analysisPlot <- renderPlot(plot(nmfTest))
        updateTabsetPanel(session, "mainTabs", selected = "analysis output")
        
    })
    
    # select which analysis to run and run it!! :) 
    observeEvent(input$runAnalysis, {
        showNotification("Running analysis. This will take a minute.")
        
        if(input$analysisChoice == "nmf"){
            print("running NMF pipeline")
            source("~/ODIS2/inst/shiny/shinyNMF.R") # issue
            globals$analysisres <- NMFforShiny(globals$dataset_list[[input$datasetRadio]], input$nmfRank)
            
            # issue -- may no longer be needed by use of sliders
            if(class(globals$analysisres) == "character"){
                showNotification(globals$analysisres)
            }
            else{
                # issue -- determine how to format res
                output$analysisPlot <- renderPlot(consensusmap(globals$analysisres))
                updateTabsetPanel(session, "mainTabs", selected = "analysis output")
                output$GSEAbutton <- renderUI({
                    tagList(checkboxGroupInput("pathwaySelection", "Select pathways for GSEA",
                                               choices = c("C1", "C2", "C3","C4","C5","C6","C7","C8", "custom (this is a placeholder and selecting it breaks the code)"),
                                               selected = NULL),
                            actionButton("runGSEA", "run GSEA"))
                })
            }
            print("NMF completed")
        } 
        # End of NMF analysis
        else{
            showNotification("not yet implemented")
            # issue -- To be completed
        }
    })
    
    # -----------------------------
    # --------- GSEA --------------
    # -----------------------------
    observeEvent(input$runGSEA, {
        showNotification("running GSEA... functioning!? :0") # :(
        print(input$pathwaySelection)
        
        source("~/ODIS2/inst/shiny/GSEAfunctions.R") # issue
        species <- switch (input$datasetRadio,
                           "LLB" = "Rattus norvegicus",
                           "snyder" = "Mus musculus",
                           "file" = "Mus musculus" # issue -- need to create some sort of input widget
        )
        print(species)
        genesets <- prepareGenesets(species, input$pathwaySelection) # issue -- allow user to select genesets. For now just testing on C8...
        #View(res)
        gene_list <- sort(globals$analysisres@fit@W[,1], decreasing = TRUE) # issue -- place holder to test with NMF, need to determine res structure
        #print(gene_list)
        globals$GSEAres <- ODISGSEA_helper(gene_list = gene_list, gmtList = genesets, pval = 1) #maintains listed structure for compatibility with ODISHEATMAP, but should only ever have one list item
        #View(GSEAres)
        
        output$genesetRadios <- renderUI({
            radioButtons("viewGeneset", "Geneset to view",
                         choiceNames = names(globals$GSEAres[[1]]),
                         choiceValues = names(globals$GSEAres[[1]]),
                         select = names(globals$GSEAres[[1]])[1]
            )
        })
        print("generating table")
        # issue "Warning: Error in [[: attempt to select less than one element in
        # get1index" I believe this may due to input$viewGeneset not existing then
        # datatable is generated, then it does not continue to show the error?
        output$gseaTable <- DT::renderDataTable(globals$GSEAres[[1]][[input$viewGeneset]][["Results"]],
                                                options = list(scrollX = TRUE))
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
