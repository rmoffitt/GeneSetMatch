# converts affymetrix ID to symbol and returns LLB expression matrix

prepareLLB <- function(){
  require(annotate)
  require(rat2302.db) #affymetrix chip for liverlungbrain
  LLB <- readRDS("~/ODIS2/data/LiverLungBrain.rds")
  
  ex <- as.data.frame(LLB$ex)
  
  genes <- rownames(ex)
  ex$PROBEID <- rownames(ex)
  convert <- AnnotationDbi::select(rat2302.db, genes, c("SYMBOL","ENTREZID", "ENSEMBL")) 
  convert <- convert[match(unique(convert$SYMBOL), convert$SYMBOL),]
  ex <- merge(ex, convert, "PROBEID")
  ex <- ex[which(!is.na(ex$ENTREZID)),]
  #rownames(W) <- W$ENTREZID
  rownames(ex) <- ex$ENTREZID
  ex <- dplyr::select(ex,!c("SYMBOL","ENTREZID", "ENSEMBL", "PROBEID"))
  
  #detach(package:annotate,unload=TRUE)
  #detach(package:rat2302.db,unload=TRUE)
  return(ex)
}

prepareSnyder <- function(){
  load("~/ODIS2/data/full_snyder.RData")
  SnyderMouseGenes <- full_Snyder$featInfo$ENSEMBL
  source("~/ODIS2/R/MouseENS2MouseENT.R")
  library(biomaRt)
  tmp_snyder <- MouseENS2MouseENT(SnyderMouseGenes)
  
  full_Snyder$featInfo <- as.data.frame(full_Snyder$featInfo)
  full_Snyder$featInfo <- merge(full_Snyder$featInfo, tmp_snyder)
  
  ex <- full_Snyder$ex
  ex$ENSEMBL <- row.names(ex)
  ex <- merge(ex, tmp_snyder, by = "ENSEMBL", all = FALSE)
  
  to.keep <- which(!is.na(ex$ENTREZID))
  ex <- ex[to.keep,]
  ex = ex[!duplicated(ex$ENTREZID),]
  row.names(ex) <- ex$ENTREZID
  ex <- dplyr::select(ex, !c("ENSEMBL", "ENTREZID", "SYMBOL"))
  return(ex)
}