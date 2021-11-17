
# facilitates feeding analysisres into ODIS.HEATMAP
# sends GSEAres with original matrix for each group/rank/etc.
# ensures symbol is included, or should this be earlier in the app?
shinyHeatmapHelper <- function(analysisres, dataex){
  #reformatting, only need clusterlist
  GSEA_res <- analysisres$clusterlist
  
  # rename genelist to be origmatrix
  for(i in length(GSEA_res)){
    GSEA_res[[i]]$orig_matrix <- GSEA_res[[i]]$genelist
    GSEA_res[[i]]<-within(GSEA_res[[i]], rm(genelist))
  }
  
  source("~/ODIS2/R/ODIS.Heatmaps.R")
  ODIS.Heatmaps(GSEA_res)
  
  return()
}