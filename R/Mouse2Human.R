#' Homology Conversion Function converts mouse ENSEMBLE ID's to human ENTREZ ID'd + HGNC and MGI's
#' @export
#' @import biomaRt
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import data.table
#' @param MouseGenes is a vector of mouse ENSEMBLE ID's
#' @return nothing
#' @examples 
#' readcounts <- data$ex
#' MouseGenes <- data$featInfo$ENSEMBL
#' humanGenes <- Mouse2Human(MouseGenes)

Mouse2Human <- function(MouseGenes){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "uswest.ensembl.org")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "uswest.ensembl.org")
  genesMousetoHuman = getLDS(attributes = c("ensembl_gene_id",
                                            "mgi_symbol"), 
                             filters = "ensembl_gene_id", 
                             values = MouseGenes , 
                             mart = mouse, 
                             attributesL = c("ensembl_gene_id",
                                             "hgnc_symbol",
                                             "entrezgene_id"), 
                             martL = human, 
                             uniqueRows = TRUE)
  colnames(genesMousetoHuman) <- c("ENSEMBL", "MGI", "Human.Gene_ID", "HGNC", "ENTREZ")
  return(genesMousetoHuman) 
}