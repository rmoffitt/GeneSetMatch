#HEATMAPS NMF
#gotta find a way to display all represented pathwars in one gmt file and all leading edges from V1, V2, V3 in one plot. Use three colors. 
# prpobably gotta restucture output withing a new heatmpa function

#




ODIS.NMF.Heatmaps <- function(gsea_results){
  
  #my_palette <- colorRampPalette(c("blue", "lightgrey", "red"))(n = 299) #MAYBE NOT YET?
  plots <- list()
  
  #this loop extracts top 30genes from each gmt.file and their top 30 associated pathways (16-39)
  for(i in names(gsea_results[[1]])){#FOR EACH LEVEL 1, PULL OBJECT Vx
    thisAnalysis <- gsea_results[[i]]
    agg_leading_edge = list() #ALL LEADCING EDGE GENES OF TOP 30 PATHWAYS
    leading_sets = character() #TOP 30 PATHWAYS FROM EACH GMT FILE
    for(k in names(gsea_results)){
      
      #make matrix and then glue smaller matrices together
      theMatrix <- matrix(data = 0,
                          nrow = length(leading_sets),
                          ncol = length(c(agg_leading_edge$V1$ENTREZID, agg_leading_edge$V2$ENTREZID, agg_leading_edge$V3$ENTREZID)))
      #I have top 30 enriched pathways
      orig_matrix = gsea_results[[k]]$orig_matrix #pull out and separate fpr appropriate k every time
      theseResults <- gsea_results[[k]][[i]]$Results[1:30,] #already sorted by normalized enrichment score
      
      plots[[i]] <- list()
      
      # Take all the leadsing edge genes from theseResults
      theEdge = unique(unlist(theseResults$leadingEdge))
      # returns the row index of theEdge vs orig_matrix, which is sorted by log2FoldChange (decreasing)
      theTopGenes = theEdge[order(match(theEdge, orig_matrix$ENTREZID))][1:30]
      agg_leading_edge[[k]] = orig_matrix[orig_matrix$ENTREZID %in% theTopGenes,c("log2FoldChange","ENTREZID")]
      print(length(theTopGenes))
  
      # agg_leading_edge = c(agg_leading_edge,
      #                      theTopGenes)
      # 
      # print(length(agg_leading_edge))
      
      leading_sets = c(leading_sets, theseResults$pathway)
    }
      #pull out orix_matrix log2fold change associated with all 90 (15-40)
      
      #do this better
      
      rownames(theMatrix) <- leading_sets
      colnames(theMatrix) <- c(agg_leading_edge$V1$ENTREZID, agg_leading_edge$V2$ENTREZID, agg_leading_edge$V3$ENTREZID)
      
      #fill TheMatrix with logfoldchange values
      
      plots[[i]][[k]] <- list()
      
      plots[[i]][[k]] <- heatmap.2(theMatrix)
      
      dev.off()
  }
}


