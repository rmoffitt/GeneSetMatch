#' This function converts Affymetrix Probe IDs into ENTREZ, SYMBOL and ENSEMBL homology. Accepts and returns a layered data.frame
#' @export
#' @import annotate
#' @import rat2302.db
#' @import base
#' @import dplyr
#' @param LLB is a data.frame of rat Affymetrix ID'sT 
#' @return a data.frame with corresponding ENTREZ, SYMBOL, ENSEMBL ID's 
#' @examples 
#' prepareLLB(LLB)

prepareLLB <- function(){
  require(annotate)
  require(rat2302.db) #affymetrix chip for liverlungbrain
  LLB <- readRDS("./LiverLungBrain.rds")
  
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
  condition <- c(rep("liver", 3), rep("brain", 3), rep("lung", 3), rep("mixed", 33))
  
  LLB <- list(ex = ex, condition= condition)
  LLB$species <- "Rattus norvegicus"
  
  #detach(package:annotate,unload=TRUE)
  #detach(package:rat2302.db,unload=TRUE)
  saveRDS(LLB, file = "~/GeneSetMatch/inst/analysis/LLB/LLB Data/LLB Data.Rds")
}
