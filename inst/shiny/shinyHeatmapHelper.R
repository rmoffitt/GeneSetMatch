
# facilitates feeding analysisres into ODIS.HEATMAP
# sends GSEAres with original matrix for each group/rank/etc.
# ensures symbol is included, or should this be earlier in the app?
shinyHeatmapHelper <- function(analysisres, species){
  #reformatting, only need clusterlist
  GSEA_res <- analysisres$clusterlist
  
  # rename genelist to be origmatrix
  for(i in 1:length(GSEA_res)){
    source("~/ODIS2/inst/shiny/getEntrezIDandSymbol.R") # :(
    print(species)
    orig_matrix <- genelistWithEntrezandSymbol(GSEA_res[[i]]$orig_matrix, species)
    
    # need to move orig_matrix to end of list based on current design of ODIS.Heatmaps
    GSEA_res[[i]]$orig_matrix <- NULL
    GSEA_res[[i]]$orig_matrix <- orig_matrix
    #View(GSEA_res)
    
    # issue -- renames column name to log2FoldChange so that it matches current ODIS HEATMAP expectations
    names(GSEA_res[[i]]$orig_matrix) <- c("ENTREZID", "SYMBOL", "log2FoldChange" )
  }
  
  source("~/ODIS2/R/ODIS.Heatmaps.R")
  library(gplots) # force libraries to load while ODIS build not working
  library(pdist)
  library(MatrixGenerics)
  plots <- ODIS.Heatmaps(GSEA_res)
  #View(plots)
  return(plots)
}