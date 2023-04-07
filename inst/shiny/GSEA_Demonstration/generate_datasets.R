 # a pretty heatmap for the powerpoint
library(gplots)
LiverLungBrain <- readRDS("~/ODIS2/data/LiverLungBrain.rds")
W <- LiverLungBrain$trueW
W <- data.matrix(W)


heatmap.2(W, scale = "column")
heatmap.2(W, scale = "row", trace = "none")

# ------------------------------
# get GSEA results to use in demo shiny

analysisres <- list()
genesets <- prepareGenesets(LLB$species, c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"))


ODISGSEA_helper <- function(analysisres, gmtList, pval){
  GSEAres <- list()
  for (cluster in 1:length(analysisres$clusterlist)) {
    for(pathway in gmtList){
      #prepare stat as named list
      
      print(pathway)
      source("~/ODIS2/R/ODISGSEA.R") # issue
      genelist <- analysisres$clusterlist[[cluster]]
      
      res <- ODISGSEA(gene_list = genelist, theGoFile = pathway, pval = pval)
      GSEAres[[names(analysisres$clusterlist)[[cluster]]]][[pathway]] <- res
    }
  }
  
  return(GSEAres)
}

for(scorename in c("per", "pm", "fc", "diff")){
  analysisres$clusterlist <- scoreNMF(LLB$genelist, scorename)
  
  LLB$GSEA[[scorename]] <- ODISGSEA_helper(analysisres = analysisres, gmtList = genesets, pval = 1)
  
  #LLB$GSEA[[scorename]] <- "place holder"
}
#saveRDS(LLB, file = "~/ODIS2/inst/shiny/GSEA_Demonstration/LLB.Rds")


# ---------------------------

#create orig_matrix
for(scorename in c("per", "pm", "fc", "diff")){
  analysisres$clusterlist <- scoreNMF(LLB$genelist, scorename)
  
  LLB$GSEA[[scorename]]$V1$orig_matrix <- analysisres$clusterlist$V1
  LLB$GSEA[[scorename]]$V2$orig_matrix <- analysisres$clusterlist$V2
  LLB$GSEA[[scorename]]$V3$orig_matrix <- analysisres$clusterlist$V3
  

}


# ------------------

shinyHeatmapHelper <- function(GSEA_res, species){
  plots <- list()
  
  for(scorename in names(GSEA_res)){
    # rename genelist to be origmatrix
    for(i in 1:length(GSEA_res[[scorename]])){
      print(species)
      orig_matrix <- genelistWithEntrezandSymbol(GSEA_res[[scorename]][[i]]$orig_matrix, species)
      
      # need to move orig_matrix to end of list based on current design of ODIS.Heatmaps
      GSEA_res[[scorename]][[i]]$orig_matrix <- NULL
      GSEA_res[[scorename]][[i]]$orig_matrix <- orig_matrix
      #View(GSEA_res)
      
      # issue -- renames column name to log2FoldChange so that it matches current ODIS HEATMAP expectations
      names(GSEA_res[[scorename]][[i]]$orig_matrix) <- c("ENTREZID", "SYMBOL", "log2FoldChange" )
    }
    
    
    source("~/ODIS2/R/ODIS.Heatmaps.R")
    library(gplots) # force libraries to load while ODIS build not working
    library(pdist)
    library(MatrixGenerics)
    plots[[scorename]] <- ODIS.Heatmaps(GSEA_res[[scorename]])
  }
  return(plots)
}


plots <- shinyHeatmapHelper(LLB$GSEA, LLB$species)
LLB$heatmapPlots <- plots
#saveRDS(LLB, file = "~/ODIS2/inst/shiny/GSEA_Demonstration/LLB.Rds")



# ===================================================================================================
# ===================================================================================================

#orig matrix for tm
require(org.Mm.eg.db)
db = org.Mm.eg.db

for(score in names(GSEA)){
  for(group in names(GSEA[[score]])){
    origMatrix <- GSEA[[score]][[group]]
    
    df <- clusterProfiler::bitr(origMatrix$gene,
                                #df <- clusterProfiler::bitr(names(origMatrix),
                                fromType = "SYMBOL", toType = "ENTREZID",
                                OrgDb = db, drop = FALSE)
    
    df$stat <- origMatrix$avg_log2FC
    
    GSEA[[score]][[group]] <- list(orig_matrix = df)
  }
}
TabulaMuris$GSEA <- GSEA

# -------------------------------------------
saveRDS(TabulaMuris_LLB, file = "~/ODIS2/inst/shiny/GSEA_Demonstration/TabulaMuris_LLB.Rds")

ODISGSEA_helper <- function(analysisres, gmtList, pval){
  GSEAres <- analysisres
  for (cluster in names(analysisres)) {
    print(cluster)
    for(pathway in gmtList){
      #prepare stat as named list
      
      print(pathway)
      source("~/ODIS2/R/ODISGSEA.R") # issue
      tmp <- analysisres[[cluster]]
      tmp <- tmp$orig_matrix[is.na(tmp$orig_matrix$ENTREZID) == F,]
      tmp <- tmp[order(tmp$stat, decreasing = T),]
      genelist <- tmp$stat
      names(genelist) <- tmp$ENTREZID
      #genelist <- sort(genelist)
      
      res <- ODISGSEA(gene_list = genelist, theGoFile = pathway, pval = pval)
      GSEAres[[cluster]][[pathway]] <- res
    }
  }
  
  return(GSEAres)
}


genesets <- prepareGenesets("Mus musculus", c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"))
GSEA <- TabulaMuris_LLB$GSEA
GSEA_res <- list()

for(scorename in names(GSEA)){
  GSEA_res[[scorename]] <- ODISGSEA_helper(analysisres = GSEA[[scorename]], gmtList = genesets, pval = 1)
}
TabulaMuris_LLB$GSEA_res <- GSEA_res

# -----------------------------------------------
tm <- TabulaMuris_LLB
library(Seurat)
Idents(tm$object) <- "cell_ontology_class"

assay <- GetAssayData(tm$object)
genelist <- data.frame(matrix(nrow = nrow(assay), ncol = 0))
for(i in unique(Idents(tm$object))){
  i_Ident <- assay[,which(Idents(tm$object) == i)]
  print(ncol(i_Ident))
  genelist[,i] <- apply(i_Ident, 1, mean)
}
#colnames(genelist) <- unique(Idents(tm$object))
rownames(genelist) <- rownames(assay)

tm$ex <- genelist
saveRDS(tm, file = "~/ODIS2/inst/shiny/GSEA_Demonstration/TabulaMuris_LLB.Rds")

# -----------------------------------

# heatmaps for Tabula Muris 

shinyHeatmapHelper <- function(GSEA_res, species = "Mus musculus"){
  plots <- list()
  
  for(scorename in names(GSEA_res)){
    # rename genelist to be origmatrix
    for(i in 1:length(GSEA_res[[scorename]])){
      #print(species)
      orig_matrix <- GSEA_res[[scorename]][[i]]$orig_matrix
      
      # need to move orig_matrix to end of list based on current design of ODIS.Heatmaps
      GSEA_res[[scorename]][[i]]$orig_matrix <- NULL
      GSEA_res[[scorename]][[i]]$orig_matrix <- orig_matrix
      #View(GSEA_res)
      
      # issue -- renames column name to log2FoldChange so that it matches current ODIS HEATMAP expectations
      names(GSEA_res[[scorename]][[i]]$orig_matrix) <- c("SYMBOL", "ENTREZID", "log2FoldChange" )
    }
    
    
    source("~/ODIS2/R/ODIS.Heatmaps.R")
    library(gplots) # force libraries to load while ODIS build not working
    library(pdist)
    library(MatrixGenerics)
    plots[[scorename]] <- ODIS.Heatmaps(GSEA_res[[scorename]])
  }
  return(plots)
}


plots <- shinyHeatmapHelper(tm$GSEA_res)
tm$heatmapPlots <- plots


saveRDS(tm, file = "~/ODIS2/inst/shiny/GSEA_Demonstration/TabulaMuris_LLB.Rds")
