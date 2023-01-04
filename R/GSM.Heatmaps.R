# GENE VS GENE SET HEATMAPS
#' This function takes direct GSEA output and transforms it into hierarchically clustered gene vs gene set heatmaps
#' @export
#' @import pdist
#' @param gsea_results is the direct output from GSEA function in step immediately prior. There is no need for restructuring, just feed it directly into the function as is. 
#' @return pdf heatmap plots corresponding to each gsea_results$Results 
#' @examples 
#' snyder_outputs is GSEA result
#' ODIS.Heatmaps(snyder_outputs)


GSM.Heatmaps <- function(gsea_results){
  
  my_palette <- colorRampPalette(c("blue", "lightgrey", "red"))(n = 299)
  plots <- list()
  filelist <- list()
  
  for(i in names(gsea_results)){
    thisAnalysis <- gsea_results[[i]]
    plots[[i]] <- list()
    for(k in (1:(length(thisAnalysis) - 1))){
      theseResults <- thisAnalysis[[k]]$Results
      
      if (length(theseResults$leadingEdge) > 40){
        theseResults <- theseResults[c(order(theseResults$ES)[1:30],
                                       order(-theseResults$ES)[1:30]),]
      }
      
      
      all.genes <- unique(c(unlist(theseResults$leadingEdge)))
      
      theMatrix <- matrix(data = 0,
                          nrow = length(theseResults$leadingEdge),
                          ncol = length(all.genes))
      
      rownames(theMatrix) <- theseResults$pathway
      
      theIndexWeWant <- match(
        as.character(all.genes),
        thisAnalysis$orig_matrix$ENTREZID) #this should be thisAnalysis$`./c1.all.v7.1.entrez.gmt`$Results$leadingEdge
      colnames(theMatrix) <- thisAnalysis$orig_matrix$SYMBOL[theIndexWeWant]
      
      plots[[i]][[k]] <- list()
      
      if(length(theseResults$leadingEdge) > 1) {
        for(j in 1:length(theseResults$leadingEdge)){
          print(j)
          theseLeadingGenes <- theseResults$leadingEdge[[j]]
          if (length(theseLeadingGenes) > 0){
            
            #where are these in my heatmap?
            firstIndexWeWant <- match(as.character(theseLeadingGenes),
                                      thisAnalysis$orig_matrix$ENTREZID[theIndexWeWant])
            
            #where are these in original data frame?
            secondIndexWeWant <- match(as.character(theseLeadingGenes), 
                                       thisAnalysis$orig_matrix$ENTREZID)
            
            theMatrix[j, firstIndexWeWant] <- thisAnalysis$orig_matrix$log2FoldChange[secondIndexWeWant] 
            }
        }
        
        print("dim of Matrix before sorting keepers")
        print(dim(theMatrix))
        
        #Here is my full matrix
        if(ncol(theMatrix) > 80){
          print("Too many columns!")
          movers = as.numeric(apply(abs(theMatrix),1,which.max)) #returns col index of largest moving gene PER ROW
          overall = abs(colSums(theMatrix)) #takes colSums
          print(range(overall))
          overall = sort(overall, decreasing = T)[1:80] #sorts big to small and takes top 80
          keepers = unique(c(which(colnames(theMatrix) %in% names(overall)), movers)) # combines top movers from each set ('movers') with consensus genes ('overall'), unique() removes overlap
          theMatrix = theMatrix[, keepers] #only keeps columns that are present in 'overall'
          
          # I've now ensured all rows have at least one non-zero, so the below is now irrelevant??
          ind = which(rowSds(theMatrix) == 0) # Finds which rows are now all 0s (sd == 0)
          print(paste0("Num Rows with zero SD: ", length(ind)))
          #theMatrix = theMatrix[ind,] # Keeps only rows where rowSD is not equal to 0, allows for clustering
        }
        print("dim of Matrix after sorting")
        print(dim(theMatrix))
        
        print(paste0(i, "_", substr(names(thisAnalysis)[k],
                                    5, stop = 999), ".pdf"))
        filelist[[i]][[k]] <- paste0(i, "_", substr(names(thisAnalysis)[k],
                                                    5, stop = 999), ".pdf")
        
        pdf(paste0(i, "_", substr(names(thisAnalysis)[k],
                                  5, stop = 999), ".pdf"),
            width = ncol(theMatrix)/10 + 5,
            height= nrow(theMatrix)/10 + 5)
        
        
        
        
        plots[[i]][[k]] <- heatmap.2(theMatrix,
                                     main = paste0(i, "_", substr(names(thisAnalysis)[k], 5, stop = 999)),
                                     distfun =  function(x) {
                                       d <- (1-cor(t(x))) + as.matrix(pdist(sign(rowMeans(x))))
                                       return(as.dist(d))},
                                     density.info = "histogram",
                                     trace = "none",
                                     margins = c(3,20),
                                     cexRow = 0.7,
                                     cexCol = 0.5,
                                     col = my_palette,
                                     breaks = (-149:150-0.5)/50,
                                     dendrogram = "both",
                                     #offsetRow = -205,
                                     #lmat = rbind(c(1,2), c(3,4), c(5,6)),
                                     #labRow = theseResults$leadingEdge,
                                     lhei = c(2, nrow(theMatrix)/20 + 3),
                                     lwid = c(2, ncol(theMatrix)/10 + 3),
                                     ##lhei=c(2,4,0.2),
                                     Colv = "Rowv"
                                     
                                     
        )
        
        dev.off()
      }
    }
  }
  
  return(filelist)
}