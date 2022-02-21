# This function generates a Hierarchically clustered Multilevel Heatmap which displays three samples in one concise plot.
#' This function takes direct GSEA output and transforms it into a unique visualization
#' @export
#' @import pdist
#' @import grid
#' @import ggplot2
#' @import sp
#' @import raster
#' @import plotly
#' @import matrixStats
#' @param gsea_results is the direct output from GSEA function in step immediately prior. There is no need for restructuring, just feed it directly into the function as is. 
#' @return multi-colored heatmaps displaying 3 experiemntal conditions/samples/vectors
#' @examples 
#' gsea_results is GSEA result
#' ODIS.Multilevel.Heatmaps(gsea_results)


#TO FIX
#####figure out "go next"
#####Run through all plots, save all
##fftw library needs updating

res_nmf <- readRDS("./Snyder_NMF_V123.rds")
res_nmf_org <- readRDS("./Snyder_NMF_ORG_V123.rds")
gsea_results <- res_nmf
i <- names(gsea_results[[1]][2])


###
ODIS.Multilevel.Heatmaps <- function(gsea_results){
  
  #define distance function for clustering
  distfun =  function(x) { 
    d <- (1-cor(t(x)))
    return(as.dist(d))}
  
  #original distance fucntion
  #library(pdist)
  # distfun1 =  function(x) {
  #   d <- (1-cor(t(x))) + as.matrix(pdist(sign(rowMeans(x))))
  #   return(as.dist(d))}
  
  #initiate list of plots
  plots <- list()
  filelist <- list()
  
  #this loop glues the small matrixes together
  for(i in names(gsea_results[[1]])){ #for each experiment V, pull object Vx 
    thisAnalysis <- gsea_results[[i]]
    theBigMatrix = list() #every time I make a small matrix for each experiment (k), glue them together here
    
    #this loop extracts top 30genes from each gmt.file and their top 30 associated pathways, and puts them into a clustered matrix
    for(k in names(gsea_results)){
      orig_matrix = gsea_results[[k]]$orig_matrix #pull out and separate for appropriate k every time
      theseResults <- gsea_results[[k]][[i]]$Results
      # if (length(theseResults == 0)){
      #   next
      # }
      theseResults <- theseResults[1:min(30, nrow(theseResults)),] #Selects up to top 30 enriched pathways and their stats by normalized enrichment score
      theEdge = unique(unlist(theseResults$leadingEdge)) #Takes all the leading edge genes from theseResults
      theTopGenes = theEdge[order(match(theEdge, orig_matrix$ENTREZID))][1:min(30, length(theEdge))] #Returns the row index of theEdge vs orig_matrix, which is sorted by log2FoldChange (decreasing)
      print(length(theTopGenes)) #checks lenght, should be 30
      agg_leading_edge = orig_matrix[orig_matrix$ENTREZID %in% theTopGenes,c("log2FoldChange","ENTREZID")] #Takes the genes corresponding to the index of theTopGenes
      print(dim(agg_leading_edge)) #checks dim of agg_leading_edge, should be 30 x 2
      
      #make a small matrix of agg_leading_edge and theseResults$pathways for each k
      theMatrix <- matrix(data = 1,
                          byrow = T,
                          nrow = length(theseResults$pathway),
                          ncol = length(agg_leading_edge$ENTREZID))
      
      print(dim(theMatrix)) #checks dim of theMatrix, should be 30X30 or so
      
      rownames(theMatrix) <- theseResults$pathway #assign rownames (pathways)
      colnames(theMatrix) <- agg_leading_edge$ENTREZID #assign colnames (genes)
      
      #for each of the leading edges, clean out theMatrix
      
      for(j in 1:length(theseResults$leadingEdge)){
        tle <- theseResults$leadingEdge[[j]]
        theMatrix[j, !(colnames(theMatrix) %in% tle)] <- 0
      }
      image(t(theMatrix))
      
      #take standard deviation of the rows and columns
      theMatrix_row_std = rowSds(theMatrix) # create a vector with all rowSds, keeps the order of rows in theMatrix
      theMatrix_col_std = colSds(theMatrix)
      
      #remove all rows and cols that have less stdev that 0.01 to avid empty space in plot
      theMatrix <- theMatrix[!(theMatrix_row_std < 0.01), !(theMatrix_col_std < 0.01)]
      
      #somewhox .x and .y is added to genenames in this merge, because genes show up in more than one V, once this is fixed 129 should work
      if(typeof(theBigMatrix)=="list"){
        #first iteration
        theBigMatrix <- theMatrix
      }else{
        theBigMatrix <- merge(theBigMatrix,theMatrix,by = "row.names",all=TRUE)
        rownames(theBigMatrix) <- theBigMatrix$Row.names #adding pathway names into row names
        theBigMatrix <- theBigMatrix[,!(names(theBigMatrix) %in% "Row.names")] #removing column with pathways names
        theBigMatrix <- as.matrix(theBigMatrix)
      }
      
      print(k)
      print(typeof(theBigMatrix)) 
      print(str(theBigMatrix))
      print(dim(theBigMatrix))
    }
    
    theBigMatrix[!is.finite(theBigMatrix)] <- 0 #removing NA/NAN/Inf from matrix before plotting
    image(t(theBigMatrix))
    
    redundantX <- grep(pattern = ".x", colnames(theBigMatrix), value = T)
    redundantX <- gsub(pattern = "\\.x", replacement = "", x = redundantX)
    
    for (z in redundantX){
      indX <- grep(colnames(theBigMatrix), pattern = z) #find the col indices of .x and .y
      theBigMatrix[which(theBigMatrix[,indX[2]] == 1), #"if there is a 1 value in the ".y", make it a 1 in the ".x"
                   indX[1]] = 1
      theBigMatrix = theBigMatrix[, -indX[2]] #remove the .y column
      colnames(theBigMatrix)[indX[1]] = z #remove the .x by replacing it with original z
    }
    
    #re-order theBigMatrix ()
    rowCluster <- hclust(distfun(theBigMatrix))
    colCluster <- hclust(distfun(t(theBigMatrix)))
    
    image(t(theBigMatrix))  
    theBigMatrixOrdered <- theBigMatrix[rowCluster$order, colCluster$order] #ordered rows and columns by rowClust$order
    image(t(theBigMatrixOrdered))
    
    #This forloop determines RGB triplet for each gene, builds corresponding matrixes and plots a raster
    rbgmatrix <- matrix(data = "#FFFFFF",
                        ncol = length(colnames(theBigMatrixOrdered)),
                        nrow = length(rownames(theBigMatrixOrdered)))
    
    colnames(rbgmatrix) = colnames(theBigMatrixOrdered)
    rownames(rbgmatrix) = rownames(theBigMatrixOrdered)
    
    for(setID in 1:nrow(theBigMatrixOrdered)){
      for(geneID in 1:ncol(theBigMatrixOrdered)){
        tbmo <- theBigMatrixOrdered[setID, geneID]
        print(tbmo)
        
        if (tbmo == 0){
          
          rbgmatrix[setID, geneID] <- "#000000" 
        }
        else {e <- colnames(theBigMatrixOrdered)[geneID]
        r <- res_nmf$V1$orig_matrix[e,"log2FoldChange"] 
        g <- res_nmf$V2$orig_matrix[e,"log2FoldChange"]
        b <- res_nmf$V3$orig_matrix[e,"log2FoldChange"]
        rbgmatrix[setID, geneID] <- rgb(r, g, b)
        }
      }
    }
    
    
    #make tall table of x, y, color for ggplot, use plotly 
    tallmatrix <- as.character(rbgmatrix)
    tallmatrix <- data.frame(data = tallmatrix,
                             x = rep(1:dim(rbgmatrix)[2], each = dim(rbgmatrix)[1]),
                             y = rep(1:dim(rbgmatrix)[1], times = dim(rbgmatrix)[2]),
                             stringsAsFactors = F)
    
    xlabels <- c(colnames(rbgmatrix))
    ylabels <- c(rownames(rbgmatrix))
    
    print(paste0(i, "_", substr(names(thisAnalysis)[k],
                                5, stop = 999), ".pdf"))
    filelist[[i]][[k]] <- paste0(i, "_", substr(names(thisAnalysis)[k],
                                                5, stop = 999), ".pdf")
    
    pdf(paste0(i, "_", substr(names(thisAnalysis)[k],
                              5, stop = 999), ".pdf"),
        width = ncol(rbgmatrix)/10 + 5,
        height= nrow(rbgmatrix)/10 + 5)
    
    
    plots[[i]][[k]] <- ggplot(tallmatrix) +
      geom_raster(aes(x = x, y= y, fill = data)) +
      scale_fill_identity() +
      scale_y_continuous(labels = ylabels, breaks = 1:length(ylabels)) +
      scale_x_continuous(labels = xlabels, breaks = 1:length(xlabels)) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
      ggtitle(paste0(i, "_", substr(names(thisAnalysis)[k], 1, stop = 999)))
    
    
    dev.off()
  }
  return(filelist)
}
