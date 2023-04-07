# selects and run DE function based on input from shiny app based on user selection
# issue -- sources all functions. need to build and call

# returns analysisres, which needs to have an object that can be plotted... or the plot somehow,
# plus clusterlist, which has the ordered list for each condition for GSEA.
DEforShiny <- function(ex, condition, method){
  #library(EnhancedVolcano)
  #library(DESeq2)
  if(method == "deseq"){
    source("~/ODIS2/R/runDESeq2_indFilterActive.R")
    deseq <- runDESeq2_indFilterActive(ex, condition)
    
    res <- list()
    res$object <- deseq$test
    #res <- res$summary
    # code below needs to be properly completed
    # hard coding 1 is problematic
    conditionNames <- levels(as.factor(condition))
    #res$clusterlist[[conditionNames[[1]]]]$genelist <- (1- deseq$summary$pval * deseq$summary$logFC)
    #names(res$clusterlist[[conditionNames[[1]]]]$genelist) <- deseq$summary$ID
    
    res$clusterlist[[conditionNames[[1]]]]$orig_matrix <- deseq$summary[c("logFC", "ID")]
    names(res$clusterlist[[conditionNames[[1]]]]$orig_matrix) <- c("stat", "ENTREZID")
    
    res$plot <- shiny::renderPlot(DESeq2::plotMA(res$object))
  }
  else if(method == "edger"){
    source("~/ODIS2/R/runedgeR_exact.R")
    res$object <- runedgeR_exact(ex, condition)$test
  }
  else if(method ==  "limma"){
    
  }
  else if(method == "etc"){
    # etc
  }
  
  # better volcano plot visual
  #output$analysisPlot <- renderPlot(EnhancedVolcano(globals$analysisres$test,
  #                                                  lab = rownames(globals$analysisres$test),
  #                                                  x = 'log2FoldChange',
  #                                                  y = 'pvalue'))

  return(res)
}
