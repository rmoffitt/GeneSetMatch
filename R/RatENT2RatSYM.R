#' Homology Conversion Function converts rat ENTREZID's to rat SYMBOL
#' @export
#' @import ClusterProfiler
#' @import data.table
#' @import msigdbr
#' @import org.Rn.eg.db
#' @param RatGenes is a vector of rat ENTREZID
#' @return nothing
#' @examples 
#' readcounts <- data$ex
#' RatGenes <- data$featInfo$ENTREZ
#' Symbols <- RatENT2RatSYM(RatGenes)

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