#' Convert rat ENTREZID's to rat SYMBOL
#' @export
#' @import ClusterProfiler
#' @import data.table
#' @import msigdbr
#' @import org.Rn.eg.db
#' @param MouseGenes is a vector of rat ENTREZID
#' @return nothing
#' @examples 
#' readcounts <- data$ex
#' MouseGenes <- data$featInfo$ENTREZ
#' Symbols <- MouseENT2MouseSYM(MouseGenes)

RatENT2RatSYM <- function(RatGenes){
  library(clusterProfiler)
  library(msigdbr)
  library(org.Rn.eg.db)
  Symbols <- bitr(RatGenes, fromType = "ENTREZID",
                  toType = c("SYMBOL"),
                  OrgDb = org.Rn.eg.db,
                  drop = FALSE) #if NA, keep entrez ID?
  return(Symbols[,2])
}