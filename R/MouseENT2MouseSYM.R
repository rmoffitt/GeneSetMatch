#' Convert mouse ENTREZID's to mouse SYMBOL
#' @export
#' @import ClusterProfiler
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

MouseENT2MouseSYM <- function(MouseGenes){
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
  Symbols <- bitr(MouseGenes, fromType = "ENTREZID",
                                toType = c("SYMBOL"),
                                OrgDb = org.Mm.eg.db)
  return(Symbols[,2])
}
