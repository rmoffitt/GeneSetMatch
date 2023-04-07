#' Convert mouse ENTREZID's to mouse SYMBOL
#' @export
#' @import clusterProfiler 
#' @import data.table
#' @import msigdbr
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @param MouseGenes is a vector of mouse ENTREZID
#' @return nothing
#' @examples 
#' readcounts <- data$ex
#' MouseGenes <- data$featInfo$ENTREZ
#' Symbols <- MouseENT2MouseSYM(MouseGenes)

MouseSYM2MouseENT <- function(MouseGenes){
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
  Symbols <- bitr(MouseGenes, fromType = "SYMBOL",
                  toType = c("ENTREZID"),
                  OrgDb = org.Mm.eg.db,
                  drop = FALSE)
  return(Symbols[,2])
}