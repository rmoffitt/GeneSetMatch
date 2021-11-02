#To do:
#make it such that it selects 'up to top 30' not 'top 30' for both genes and pathways inside k loop
#create colors based on K and not my_pallette
#add rownames to theBigMatrix (pathways)
#make heatmap
#line 76?



ODIS.NMF.Heatmaps <- function(gsea_results){
  
  #define distance function for clustering
  distfun =  function(x) { 
    d <- (1-cor(t(x)))
    return(as.dist(d))}
  
  #define color palette - different color for each K and range of pigment for log2foldchange
  my_palette <- colorRampPalette(c("blue", "lightgrey", "red"))(n = 299)
  
  #initiate list of plots
  plots <- list()
  
  #this loop glues the small matrixes together
  for(i in names(gsea_results[[1]])){ #for each experimnet V, pull object Vx 
    thisAnalysis <- gsea_results[[i]]
    theBigMatrix = list() #every time I make a small matrix for each experiemnt (k), glue them together here
    
    #this loop extracts top 30genes from each gmt.file and their top 30 associated pathways, and puths them into a clustered matrix
    for(k in names(gsea_results)){
      orig_matrix = gsea_results[[k]]$orig_matrix #pull out and separate for appropriate k every time
      theseResults <- gsea_results[[k]][[i]]$Results[1:25,] #Selects top 30 enriched pathways and their statsby normalized enrichment score
      theEdge = unique(unlist(theseResults$leadingEdge)) #Takes all the leading edge genes from theseResults
      theTopGenes = theEdge[order(match(theEdge, orig_matrix$ENTREZID))][1:200] #Returns the row index of theEdge vs orig_matrix, which is sorted by log2FoldChange (decreasing)
      print(length(theTopGenes)) #checks lenght, should be 30
      agg_leading_edge = orig_matrix[orig_matrix$ENTREZID %in% theTopGenes,c("log2FoldChange","ENTREZID")] #Takes the genes corresponding to the index of theTopGenes
      print(dim(agg_leading_edge)) #checks dim of agg_leading_edge, should be 30 x 2
      
      #make a small matrix of agg_leading_edge and theseResults$pathways for each k
      theMatrix <- matrix(data = agg_leading_edge$log2FoldChange,
                          byrow = T,
                          nrow = length(theseResults$pathway),
                          ncol = length(agg_leading_edge$ENTREZID))

      print(dim(theMatrix)) #checks dim of theMatrix, should be 30X30 or so
      
      rownames(theMatrix) <- theseResults$pathway #assign rownames (pathways)
      colnames(theMatrix) <- agg_leading_edge$ENTREZID #assign colnames (genes)
        
      #for each of the leading edges, clean out theMatrix
      for(j in 1:length(theseResults$leadingEdge)){
        tle <- theseResults$leadingEdge[[j]]
        theMatrix[j, !(colnames(theMatrix) %in% tle)] <- -0.01
      }
      image(t(theMatrix))
      #Cluster theMatrix by distfun (defined at beginning)
      rowCluster <- hclust(distfun(theMatrix))
      colCluster <- hclust(distfun(t(theMatrix)))
      
      if(typeof(theBigMatrix)=="list"){
        #first iteration
        theBigMatrix <- theMatrix
      }else{
        theBigMatrix <- merge(theBigMatrix,theMatrix,by = "row.names",all=TRUE)
        # get rid of weird columns, then turn it back into a matrix
        theBigMatrix <- theBigMatrix[,!(names(theBigMatrix) %in% "Row.names")]
        #add rownames(pathways) to TheBigMatrix
        theBigMatrix <- as.matrix(theBigMatrix)
      }
      print(k)
      print(typeof(theBigMatrix)) 
      print(str(theBigMatrix))
      print(dim(theBigMatrix))
    }
image(t(theBigMatrix))    
image(theBigMatrix)
theBigMatrix[!is.finite(theBigMatrix)] <- 0 #removing NA/NAN/Inf from matrix befroe plotting

# heatmap it, it should already be clustered
      
      plots[[i]][[k]] <- list()
      
      pdf(paste0(i, "_", substr(names(thisAnalysis)[k],
                                5, stop = 999), ".pdf"),
          width = ncol(theBigMatrix)/10 + 5,
          height= nrow(theBigMatrix)/10 + 5)
      
      plots[[i]][[k]] <- heatmap.2(theBigMatrix,
                                   main = paste0(i, "_", substr(names(thisAnalysis)[k], 5, stop = 999)),
                                   hclustfun = distfun,
                                   #distfun = distfun,
                                   density.info = "histogram",
                                   dendrogram = "none",
                                   trace = "none",
                                   margins = c(3,20),
                                   cexRow = 0.7,
                                   cexCol = 0.5,
                                   col = my_palette, #color by thisAnalysis[k]
                                   #breaks = (-149:150-0.5)/50,
                                   #dendrogram = "both",
                                   #offsetRow = -205,
                                   #lmat = rbind(c(1,2), c(3,4), c(5,6)),
                                   #labRow = theseResults$leadingEdge,
                                   lhei = c(2, nrow(theBigMatrix)/20 + 3),
                                   lwid = c(2, ncol(theBigMatrix)/10 + 3),
                                   ##lhei=c(2,4,0.2),
                                   Colv = "Rowv"
                                   
                                   
      )
      
      
      dev.off()
  }
}


