#
# Stripped down version of app to demonstrate already run functions for Capstone presentation
# 

library(shiny)
library(NMF)
library(Seurat)
options(shiny.maxRequestSize = 30*1024^2) #30MB
source("~/ODIS2/inst/shiny/GSEA_Demonstration/helper.R")

# currently user can change selection throughout the process, potentially creating strange results... maybe save selections or prevent user from changing selection after each step.

# -------------------------------------------------------------
# ---------------------------- UI -----------------------------
# -------------------------------------------------------------
ui <- fluidPage(
    # Application title
    titlePanel("Ranking for GSEA"),
    
    # --------------
    # -- SIDE BAR --
    # --------------
    column(3,
           wellPanel(
               h4("In this demonstartion you can view the liver/lung/brain bulk RNAseq or Tabula Muris scRNAseq"),
               # dataset UI
               radioButtons("datasetRadio", "select dataset",
                            choices = list("Liver/Lung/Brain" = "LLB", "Tabula Muris" = "tm"),
                            selected = "LLB"),
           ),
           
        #=====================================================================================
           # analysis UI
           wellPanel(
               # option will be dictated by which dataset has been selected
               h3("Select Analysis"),
               radioButtons("analysisChoice", "Analysis",
                            choices = list("NMF" = "nmf", "FindAllMarkers" = "fam"),
                            selected = "nmf"),
               
               p("-----------------------", align = "center"),
               # Create options based on method selection
               conditionalPanel("input.analysisChoice == 'nmf'",
                                p("In this demo, rank is set to 3."),
                                sliderInput("nmfRank", "rank to use:",
                                            min = 3, max = 3,
                                            value = 3,
                                            ticks = F
                                )
               ),
               radioButtons("analysisOptions", "Select list ranking option",
                            choices = list("none" = "none"),
                            selected = "none")
           ),
           #============================================================================================
           # gsea options UI
           wellPanel(
               uiOutput("clusterRadioUI"),
               uiOutput("genesetRadios")
            )
    ),
    
    
    # ---------------------
    # -- MAIN PANEL TABS --
    # ---------------------
    column(9,
           tabsetPanel(id = "mainTabs",
                       type = "tabs",
                       tabPanel("Preview raw data", p("Average expression"), 
                                DT::dataTableOutput("rawData")),
                       # issue -- set plot output size to prevent weird stretching
                       tabPanel("analysis output", plotOutput("analysisPlot"), 
                                uiOutput("clusterRadio"), 
                                DT::dataTableOutput("orderedData")
                                ),
                       tabPanel("GSEA output",
                                DT::dataTableOutput("gseaTable")
                                ),
                       tabPanel("ODIS Heatmap", uiOutput("heatmapPDF")),
                       selected = "Preview raw data"
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
    globals <- reactiveValues(dataset_list = list(LLB = readRDS("~/ODIS2/inst/shiny/GSEA_Demonstration/LLB.Rds"),
                                                  tm = readRDS("~/ODIS2/inst/shiny/GSEA_Demonstration/TabulaMuris_LLB.Rds"),
                                                  file = NULL), 
                              analysisres = NULL,
                              #GSEAres = NULL, # no longer used, just stores in analysisres
                              heatmapPlots = NULL
    )
    
    # -------------
    # -- DATASET --
    # -------------
    
    # console debug helper, prints selected dataset into console when selection changes.
    observeEvent(input$datasetRadio, {
        #print(input$datasetRadio)
        if(input$datasetRadio == "LLB"){
            updateRadioButtons(session, "analysisChoice",
                               choiceNames = c("NMF"),
                               choiceValues = c("nmf"),
                               selected = "nmf")
        } else if(input$datasetRadio == "tm"){
            updateRadioButtons(session, "analysisChoice",
                               choiceNames = c("FindAllMarkers"),
                               choiceValues = c("fam"),
                               selected = "fam")
        }
    })

    
    #set condition for input dataset
    observeEvent(input$conditionLab, {
        globals$dataset_list$file$condition <- unlist(strsplit(input$conditionLab, ","))
    })
    #set species for input dataset
    observeEvent(input$fileSpecies, {
        globals$dataset_list$file$species <- input$fileSpecies
    })
    
    # creates table of dataset that is displayed in raw data tab
    output$rawData <- DT::renderDataTable(
        globals$dataset_list[[input$datasetRadio]]$ex,
        options = list(scrollX = TRUE)
    )
    
    observeEvent(input$analysisChoice, {
        #print(input$analysisChoice)
        if(input$analysisChoice == "nmf"){
            updateRadioButtons(session, "analysisOptions",
                               choiceNames = c("percent expression", "transformed pattern markers", "log2fc", "diff of cols"),
                               choiceValues = c("per", "pm", "fc", "diff"),
                               selected = "per")
        } else if(input$datasetRadio == "tm"){
            updateRadioButtons(session, "analysisOptions",
                               choiceNames = c("Deafult FindAllMarkers", "Default fill 0", "Max FindAllMarkers", "Max fill 0"),
                               choiceValues = c("default", "default0", "max", "max0"),
                               selected = "default")
        }
    })
    
    # --------------
    # -- ANALYSIS --
    # --------------
    # show plot
    observeEvent(input$datasetRadio, {
        if(input$datasetRadio == "LLB"){
            output$analysisPlot <- renderPlot(consensusmap(globals$dataset_list[[input$datasetRadio]]$analysis))
        }
        else if(input$datasetRadio == "tm"){
            library(Seurat)
            output$analysisPlot <- renderPlot(DimPlot(globals$dataset_list[[input$datasetRadio]]$object, group.by = "cell_ontology_class", reduction = "umap"))
        }
    })
    #show ordered genelist
    observeEvent(input$analysisOptions, {
        globals$dataset_list[[input$datasetRadio]]$ordered <- rankGenelist(globals$dataset_list[[input$datasetRadio]]$genelist, input$analysisChoice, input$analysisOptions)
        
        #print(names(globals$dataset_list[[input$datasetRadio]]$ordered))
        
        output$orderedData <- DT::renderDataTable(
            createDataTable(globals$dataset_list[[input$datasetRadio]]$ordered),
            options = list(scrollX = TRUE)
        )
    })
    
    # ----------
    # -- GSEA --
    # ----------
    
    #output$GSEA <- renderUI({
    #    tagList(wellPanel(
    #        h3("GSEA"),
    #        checkboxGroupInput("pathwaySelection", "Select pathways for GSEA",
    #                           choices = c("C1", "C2", "C3","C4","C5","C6","C7","C8"), #"custom (this is a placeholder and selecting it breaks the code)"),
    #                           selected = NULL),
    # fileInput for custom pathways?
    # textInput for custom pathways?
    #        actionButton("runGSEA", "run GSEA"))
    #    )
    #})
    
    
    output$genesetRadios <- renderUI({
        clusterPathwayNames <- names(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]][[input$clusterRadio]])
        clusterPathwayNames <- clusterPathwayNames[-which(clusterPathwayNames == "orig_matrix")] # remove orig_matrix
        #print(clusterPathwayNames)
        
        tagList( 
            radioButtons("viewGeneset", "Geneset to view",
                         choiceNames = clusterPathwayNames,
                         choiceValues = clusterPathwayNames,
                         select = clusterPathwayNames[1]
            )
            #,p("Heatmap will be created for currently selected cluster and geneset category")
            #,actionButton("makeHeatmap", "generate heat map")
        )
    })
    
    output$clusterRadioUI <- renderUI({
        selectInput("clusterRadio", "Cluster to view GSEA results for",
                    choices = names(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]]), 
                    selected = names(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]])[[1]])
    })
    
    
    output$gseaTable <- DT::renderDataTable(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]][[input$clusterRadio]][[input$viewGeneset]][["Results"]],
                                            options = list(scrollX = TRUE))

    
    
    
    # -------------
    # -- Heatmap --
    # -------------
    output$heatmapPDF <- renderUI({
        #print("HeatmapPDF")
        #print(names(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]][[input$clusterRadio]]))
        if(input$datasetRadio == "tm"){
        i <- which(names(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]][[input$clusterRadio]]) == input$viewGeneset) -1
        }else{
        i <- which(names(globals$dataset_list[[input$datasetRadio]]$GSEA[[input$analysisOptions]][[input$clusterRadio]]) == input$viewGeneset)
        }
        tagList(
            # issue... getwd() doesn't seem to work, so not sure how to determine where file has been saved to.
            # this bandaid works for now though...!
            tags$iframe(style="height:600px; width:100%", src=paste0("http://104.41.140.14:8787/files/ODIS2/inst/shiny/GSEA_Demonstration/heatmaps/", input$analysisOptions, "/", globals$dataset_list[[input$datasetRadio]]$heatmapPlots[[input$analysisOptions]][[input$clusterRadio]][[i]]))
        )
        })
}

# Run the application 
shinyApp(ui = ui, server = server)