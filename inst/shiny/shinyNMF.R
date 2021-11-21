library(NMF)
library(matrixStats)

# function designed to return NMF output compatible with downstream GSEA analysis
# Will return the ordered W matrix, or string describing error encountered
NMFforShiny <- function(dataset, rank = 3, scoreMethod, nrun = 30){
  # issue -- check not implemented
  #check that dataset contains only non-negative numbers
  print("removing rows where all values are 0")
  tmp <- dataset[(rowSums(dataset) == 0),]
  print(paste(nrow(tmp), "rows removed"))
  dataset <- dataset[(rowSums(dataset) != 0),]
  #print(class(dataset))
  #print(ncol(dataset))
  
  #check that rows are named? or check this somewhere else?
  #print(row.names(dataset))
  
  #check that rank is < ncol and > 1
  # these are no longer possible with the slider options.
  if(rank > ncol(dataset)){
    return("NMF not completed: rank must be less than column number.")
  }
  else if(rank <= 1){
    return("NMF not completed: rank must be greater than 1.")
  }
  
  
  #run NMF
  nmfresult <- nmf(dataset, rank, .options = "p2")
  W <- nmfresult@fit@W
  
  # methods to order W matrix
  #"percent expression" = "per", "pattern markers" = "pm", "log2fc" = "fc", "diff of cols" = "diff"
  print("ordering by W value")
  source("~/ODIS2/inst/shiny/nmf_scoring_methods.R")
  if(scoreMethod == "per"){
    wscored <- percentScore(W)
  }
  else if(scoreMethod == "pm"){
    wscored <- patternMarker(W)
  }
  else if(scoreMethod == "fc"){
    wscored <- fc(W)
  }
  else if(scoreMethod == "diff"){
    wscored <- diffAvg(W)
  }

  
  analysisres <- list()
  #analysisres$clusterlist <- vector(mode = "list", length = rank)
  for(i in 1:rank) {
    origMatrix <- data.frame(stat = wscored[, i], ENTREZID = row.names(W))
    analysisres$clusterlist[[paste0("V",i)]]$orig_matrix <- origMatrix
    #analysisres$clusterlist[[paste0("V",i)]]$genelist <- wscored[ , i]
  }
  #names(analysisres) <- colnames(nmfresult)
  analysisres$object <- nmfresult
  analysisres$plot <- renderPlot(consensusmap(analysisres$object))
  
  return(analysisres)
}

estimateRankShiny <- function(dataset, ranks){
  print("removing rows where all values are 0")
  tmp <- dataset[(rowSums(dataset) == 0),]
  print(paste(nrow(tmp), "rows removed"))
  dataset <- as.matrix(dataset[(rowSums(dataset) != 0),])
  
  showNotification("Rank test is run with 1000 variable genes to reduce run time")
  data.reduce <- dataset[order(rowVars(dataset), decreasing=TRUE),]
  data.reduce <- data.reduce[1:1000,]
  
  nmfTest <- nmf(data.reduce, ranks, nrun = 10, .options = "p2")
  return(nmfTest)
}
