# prepare tabula muris
#TabulaMuris_LLB <- readRDS("~/ODIS2/inst/shiny/GSEA_Demonstration/TabulaMuris_LLB.Rds")
#saveRDS(LLB, "~/ODIS2/inst/shiny/GSEA_Demonstration/LLB.Rds")




# -------------------------------
# Analysis helper functions

rankGenelist <- function(genelist, method, option){
  print(method)
  print(option)
  if(method == "nmf"){
    scored <- scoreNMF(genelist, option)
  }
  else if(method == "fam"){
    scored <- scoreFAM(genelist, option)
  }
  
  
  
  return(scored)
}


scoreFAM <- function(FAM, method){
  #View(FAM)
  genelist <- FAM[[method]]
  #View(genelist)
  
  clusters <- names(genelist)
  first <- clusters[[1]]
  df <- data.frame(stat = genelist[[first]]$orig_matrix$stat, gene = genelist[[first]]$orig_matrix$SYMBOL)
  colnames(df) <- c(first, "gene")
  
  clusters <- clusters[-1]
  
  for(cluster in clusters){
    tmp <- data.frame(stat = genelist[[cluster]]$orig_matrix$stat, gene = genelist[[cluster]]$orig_matrix$SYMBOL)
    tmp[[cluster]] <- tmp$stat
    tmp$stat <- NULL
    df <- merge(df, tmp, by = "gene", all = TRUE)
  }
  
  row.names(df) <- df$gene
  df$gene <- NULL
  
  return(df)
}

scoreNMF <- function(W, scoreMethod){
  #run NMF
  
  # methods to order W matrix
  #"percent expression" = "per", "pattern markers" = "pm", "log2fc" = "fc", "diff of cols" = "diff"
  #print(scoreMethod)
  source("~/ODIS2/inst/shiny/nmf_scoring_methods.R")
  if(scoreMethod == "per"){
    wscored <- percentScore(W)
  }
  else if(scoreMethod == "pm"){
    wscored <- transformPM(W)
  }
  else if(scoreMethod == "fc"){
    wscored <- fc(W)
  }
  else if(scoreMethod == "diff"){
    wscored <- diffAvg(W)
  }
  else{
    wscored <- NULL
  }
  
  if(is.null(wscored) == FALSE){
    clusters <- list()
    #analysisres$clusterlist <- vector(mode = "list", length = rank)
    for(i in 1:ncol(wscored)) {
      clusters[[paste0("V",i)]] <- wscored[,i]
    }
    
    
    wscored <- clusters
  }
  
  return(wscored)
}

# for displaying on the datatable
createDataTable <- function(ordered){
  df <- as.data.frame(ordered[1])
  for(i in 2:length(ordered)){
    df <- merge(df, as.data.frame(ordered[i]), by = "row.names")
    rownames(df) <- df$Row.names
    df <- subset(df, select = -Row.names)
  }
  
  return(df)
}
