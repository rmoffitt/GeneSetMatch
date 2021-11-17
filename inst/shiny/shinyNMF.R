library(NMF)
library(matrixStats)

# function designed to return NMF output compatible with downstream GSEA analysis
# Will return the ordered W matrix, or string describing error encountered
NMFforShiny <- function(dataset, rank = 3, scoreMethod = "log2", nrun = 30){
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
  print("nmf function")
  nmfresult <- nmf(dataset, rank, .options = "p2")
  #score W matrix
  
  #return something
  return(nmfresult)
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
