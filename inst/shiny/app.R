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

# currently user can change selection throughout the process, potentially creating strange results... maybe save selections or prevent user from changing selection after each step.

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
                         selected = "LLB"),
            
            conditionalPanel("input.datasetRadio == 'file'",
                             radioButtons("fileSpecies", "select species for input dataset",
                                          choices = list("Human" = "Homo sapiens",
                                                         "Mouse" = "Mus musculus",
                                                         "Rat" = "Rattus norvegicus")))
        ),
        #=====================================================================================
        # analysis UI
        wellPanel(
            h3("Select Analysis"),
            radioButtons("analysisChoice", "Analysis",
                         choices = list("Differential expression" = "de", "NMF" = "nmf", "FindAllMarkers" = "fam",
                                        "etc" = "null"),
                         selected = "de"),
            p("-----------------------", align = "center"),
            # Create options based on method selection
            conditionalPanel("input.analysisChoice == 'de'",
                             p("DE options go here"),
                             selectInput("DEmethod", "Method to perform DE:",
                                         choices = list("DEseq" = "deseq", "edgeR" = "edger", "etc" = "etc"),
                                         selected = "deseq"),
                             selectInput("test", "test of subcategories, does not do anything with choice:",
                                         choices = list("DEseq" = list("Cookoff", "filter", "etc"), 
                                                        "edgeR" = list("option", "option2"), 
                                                        "etc" = "etc"),
                                         selected = "deseq"),
                             conditionalPanel("input.datasetRadio == 'file'",
                                              textInput("conditionLab",
                                                        "A comma seperated list of groups/condition for each sample
                                                        in the same order as the uploaded dataset.
                                                        Ex: one,two,one,three,two,one")
                                              )
            ),
            conditionalPanel("input.analysisChoice == 'nmf'",
                             p("You can test several ranks, or simply select a rank if you already know what to use."),
                             sliderInput("rankToTest", "Ranks to test:",
                                         min = 2, max = 20, #change max to ncol? Just cap at something reasonable?
                                         value = c(3, 10),
                                         ticks = F
                             ),
                             actionButton("rankTest", "test different ranks"),
                             br(),
                             br(),
                             p("-----------------------", align = "center"),
                             sliderInput("nmfRank", "rank to use:", #issue -- ideally would not display blue on color bar, but fine for now
                                         min = 2, max = 20,
                                         value = 3,
                                         ticks = F
                             ),
                             radioButtons("nmfMethod", "select method for ordering W",
                                          choices = list("percent expression" = "per", "pattern markers" = "pm", "log2fc" = "fc", "diff of cols" = "diff"),
                                          selected = "per")
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
        # gsea options UI
        uiOutput("GSEA")
    ),

            
    # ---------------------
    # -- MAIN PANEL TABS --
    # ---------------------
    column(9,
           tabsetPanel(id = "mainTabs",
                       type = "tabs",
                       tabPanel("Preview raw data", DT::dataTableOutput("rawData")),
                       # issue -- set plot output size to prevent weird stretching
                       tabPanel("analysis output", plotOutput("analysisPlot"), textOutput("gseaPrompt")),
                       tabPanel("GSEA output", uiOutput("clusterRadioUI"), uiOutput("genesetRadios"), DT::dataTableOutput("gseaTable")),
                       tabPanel("ODIS Heatmap", uiOutput("heatmapUI"), uiOutput("heatmapPDF"))
                       #tabPanel("PDF test", tags$iframe(style="height:900px; width:75%", src="http://104.41.140.14:8787/files/ODIS2/Negative_C1.gmt.pdf"))
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
    globals <- reactiveValues(dataset_list = list(LLB = readRDS("~/ODIS2/data/LLB.Rds"), # issue -- build and set loads to ODIS::dataname
                                             snyder = readRDS("~/ODIS2/data/full_snyder_ex.Rds"),
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
        print(input$datasetRadio)
    })
    
    #gets the name of the uploaded file--> used to update radio buttons for dataset options.
    uploadedFileName <- eventReactive(input$dataFile, {
        input$dataFile$name
    })
    
    # reads in data file and updates dataset radio button options 
    observeEvent(input$dataFile, {
        updateRadioButtons(session, "datasetRadio",
                           label = paste("Upload dataset or use example dataset"),
                           choiceNames = c("Liver/Lung/Brain", "Snyder", uploadedFileName()),
                           choiceValues = c("LLB", "snyder", "file"),
                           selected = "file")
        
        #load in file
        # issue -- check what type of file, want to accept Rdata/Rds dataframes and csv
        e = new.env()
        globals$dataset_list$file$ex <- e[[load(input$dataFile$datapath, envir = e)]]
        #if(csv)
        #readcsv
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
    
    # --------------
    # -- ANALYSIS --
    # --------------
    
    #NMF rank test
    observeEvent(input$rankTest, {
        showNotification("Running nmfEstimateRank. This will take a minute.")
        ranks <- seq(input$rankToTest[1], input$rankToTest[2])
        print(ranks)
        nmfTest <- estimateRankShiny(globals$dataset_list[[input$datasetRadio]]$ex, ranks)
        print("done")
        output$analysisPlot <- renderPlot(plot(nmfTest))
        updateTabsetPanel(session, "mainTabs", selected = "analysis output")
        
    })
    
    # select which analysis to run and run it!! :) 
    observeEvent(input$runAnalysis, {
        showNotification("Running analysis. This will take a minute.")
        
        if(input$analysisChoice == "de"){
            source("~/ODIS2/inst/shiny/shinyDE.R") # issue
            print("running differential expresion pipeline")
            #print(globals$dataset_list[[input$datasetRadio]]$condition)
            globals$analysisres <- DEforShiny(globals$dataset_list[[input$datasetRadio]]$ex, globals$dataset_list[[input$datasetRadio]]$condition, input$DEmethod)
            
            print("DE done")
        }
        else if(input$analysisChoice == "nmf"){
            print("running NMF pipeline")
            source("~/ODIS2/inst/shiny/shinyNMF.R") # issue
            globals$analysisres <- NMFforShiny(globals$dataset_list[[input$datasetRadio]]$ex, input$nmfRank, input$nmfMethod)
            
            print("NMF completed")
        } 
        else{
            showNotification("not yet implemented")
            # issue -- To be completed
        }
        
        
        # change output now that analysis is completed
        #View(globals$analysisres)
        output$analysisPlot <- globals$analysisres$plot
        updateTabsetPanel(session, "mainTabs", selected = "analysis output")
        output$GSEA <- renderUI({
            tagList(wellPanel(
                h3("GSEA"),
                checkboxGroupInput("pathwaySelection", "Select pathways for GSEA",
                                   choices = c("C1", "C2", "C3","C4","C5","C6","C7","C8"), #"custom (this is a placeholder and selecting it breaks the code)"),
                                   selected = NULL),
                # fileInput for custom pathways?
                # textInput for custom pathways?
                actionButton("runGSEA", "run GSEA"))
            )
        })
        output$gseaPrompt <- renderText("Select or import pathways in the left panel and start GSEA analysis when ready.")
    })
    
    # ----------
    # -- GSEA --
    # ----------
    observeEvent(input$runGSEA, {
        showNotification("running GSEA... functioning!? :0") # :(
        print(input$pathwaySelection)
        #View(globals$analysisres)
        
        source("~/ODIS2/inst/shiny/GSEAfunctions.R") # issue
        #print(globals$dataset_list[[input$datasetRadio]]$species)
        genesets <- prepareGenesets(globals$dataset_list[[input$datasetRadio]]$species, input$pathwaySelection) # issue -- allow user to select genesets.
        globals$analysisres <- ODISGSEA_helper(analysisres = globals$analysisres, gmtList = genesets, pval = 1)
        
        
        # update UI after GSEA completes
        output$clusterRadioUI <- renderUI({
            selectInput("clusterRadio", "Cluster to view GSEA results for",
                        choices = names(globals$analysisres$clusterlist), 
                        selected = names(globals$analysisres$clusterlist)[1])
        })
        
        output$genesetRadios <- renderUI({
            clusterPathwayNames <- names(globals$analysisres$clusterlist[[input$clusterRadio]])
            clusterPathwayNames <- clusterPathwayNames[-which(clusterPathwayNames == "orig_matrix")] # remove orig_matrix
            print(clusterPathwayNames)
            
            tagList( 
            radioButtons("viewGeneset", "Geneset to view",
                         choiceNames = clusterPathwayNames,
                         #choiceValues = names(globals$analysisres$clusterlist[[input$clusterRadio]])[-1],
                         choiceValues = clusterPathwayNames,
                         select = clusterPathwayNames[1]
            ),
            actionButton("makeHeatmap", "generate heat map")
            )
        })
        
        print("generating table")
        # issue "Warning: Error in [[: attempt to select less than one element in get1index" I believe this may due to input$viewGeneset not existing when datatable is generated, then it does not continue to show the error?
        output$gseaTable <- DT::renderDataTable(globals$analysisres$clusterlist[[input$clusterRadio]][[input$viewGeneset]][["Results"]],
                                                options = list(scrollX = TRUE))
        updateTabsetPanel(session, "mainTabs", selected = "GSEA output")
    })
    
    
    
    # -------------
    # -- Heatmap --
    # -------------
    observeEvent(input$makeHeatmap, {
        showNotification("You clicked the button! :)")
        #View(globals$analysisres)
        print(globals$dataset_list[[input$datasetRadio]]$species)
        
        source("~/ODIS2/inst/shiny/shinyHeatmapHelper.R") # issue
        globals$heatmapPlots <- shinyHeatmapHelper(analysisres = globals$analysisres, species = globals$dataset_list[[input$datasetRadio]]$species)
        
        updateTabsetPanel(session, "mainTabs", selected = "ODIS Heatmap")
        # https://gist.github.com/aagarw30/d5aa49864674aaf74951
        # imbedding PDF images
        
        output$heatmapUI <- renderUI({
            #print(names(globals$heatmapPlots))
            #print(input$pathwaySelection)
            #print(length(input$pathwaySelection))
            tagList(
                radioButtons("heatmapCluster", "Select cluster to view heatmap for:",
                             choiceNames = names(globals$heatmapPlots),
                             choiceValues = names(globals$heatmapPlots),
                             selected = names(globals$heatmapPlots)[[1]]),
                radioButtons("heatmapPathway","Select pathway category to view heatmap for:",
                            choiceNames = input$pathwaySelection,
                            choiceValues = 1:length(input$pathwaySelection)),
                #            selected = 1),
                
                #output$heatmapPlot <- renderPlot(heatmap.2(globals$heatmapPlots[[input$heatmapCluster]][[input$heatmapPathway]])), #when return plots from ODIS.Heatmaps... does not seem to be a way too recreate the heatmap.2 when returned
                
                #renderText("showing", paste0("http://104.41.140.14:8787/files/ODIS2/", "V1_C8.gmt", ".pdf")),#, getwd()),
            )
        })
        
        #print(input$heatmapCluster)
        #print(input$heatmapPathway)
        
    })
    
    observeEvent(input$heatmapCluster, {
        print(input$heatmapCluster)
    })
    
    observeEvent(input$heatmapPathway, {
        print(input$heatmapPathway)
        print(globals$heatmapPlots[[input$heatmapCluster]][[as.numeric(input$heatmapPathway)]])
    })
    
    output$heatmapPDF <- renderUI(
        tagList(
            # issue... getwd() doesn't seem to work, so not sure how to determine where file has been saved to.
            # this bandaid works for now though...!
            tags$iframe(style="height:600px; width:100%", src=paste0("http://104.41.140.14:8787/files/ODIS2/inst/shiny/", globals$heatmapPlots[[input$heatmapCluster]][[as.numeric(input$heatmapPathway)]]))
            ))
}

# Run the application 
shinyApp(ui = ui, server = server)
