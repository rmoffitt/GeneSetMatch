#' Convert mouse ENSEMBLE ID's to mouse ENTREZ ID's + MGI's
#' @export
#' @import biomaRt
#' @import data.table
#' @import msigdbr
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @param MouseGenes is a vector of mouse ENSEMBLE ID's
#' @return nothing
#' @examples 
#' readcounts <- data$ex
#' mouseENS <- data$featInfo$ENSEMBL
#' mouseENT <- MouseENS2MouseENT(MouseGenes)

MouseENS2MouseENT <- function(MouseGenes){
  mouseENT = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "useast.ensembl.org")
  mouseENS = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "useast.ensembl.org")
  genesMouseENS2MouseENT = getBM(attributes = c("entrezgene_id", "ensembl_gene_id", "mgi_symbol"),
                             filters = "ensembl_gene_id", 
                             values = MouseGenes, 
                             mart = mouseENS,
                             uniqueRows = TRUE)
  colnames(genesMouseENS2MouseENT) <- c("ENTREZID", "ENSEMBL", "SYMBOL")
  return(genesMouseENS2MouseENT) 
}