# This function generates a Hierarchically clustered Multilevel Heatmaps which displays three samples in one concise plot.
#' This function takes direct GSEA output and transforms it into a unique visualization
#' @export
#' @import pdist
#' @import grid
#' @import ggplot2
#' @import sp
#' @import raster
#' @import plotly
#' @import matrixStats
#' @import DescTools
#' @param gsea_results is the direct output from GSEA function in step immediately prior. There is no need for restructuring, just feed it directly into the function as is. 
#' @return multi-colored heatmaps displaying 3 experimental conditions/samples/vectors
#' @examples 
#' gsea_results is GSEA result
#' ODIS.Multilevel.Heatmaps(gsea_results)


##DATA OF CHOICE
# res_nmf <- readRDS("./Snyder_NMF_V123.rds")
snyder_nmf_org <- readRDS("./Snyder_NMF_ORG_GSEA_NEW.rds")
#res_nmf_cell <- readRDS("./Snyder_NMF_CELL_GSEA.rds")
#res_nmf_llb_sym <- readRDS("./LLB_NMF_GSEA_SYM.rds")
#res_nmf_llb <- readRDS("./LLB_NMF_GSEA.rds")
#res_nmf_tfpm <- readRDS("./LLB_NMF_GSEA_TFPM.rds")
llb_new_order <- readRDS("./LLB_NMF_GSEA_TFPM_NEW.rds")
llb_new_pval <- readRDS("./LLB_NMF_GSEA_TFPM_NEW_0.1p.rds")

gsea_results <- llb_new_pval

i <- names(gsea_results[[3]][9])

###
ODIS.Multilevel.Heatmaps <- function(gsea_results){
  
  library(ggplot2)
  library(DescTools)
  
  #define distance function for clustering
  distfun =  function(x) { 
    d <- (1-cor(t(x)))
    return(as.dist(d))}
  
  HSV.ROW.distfun <- function(x,weight) {
    # figure the color for that row
    averageHue <- numeric(0)
    for( l in 1:(dim(x)[1])){
      thisrow <- x[l,]
      h.x <- cos(thisrow*2*pi)
      h.y <- sin(thisrow*2*pi)
      h <- complex(real = h.x , imaginary = h.y)
      h.mean <- mean(h,na.rm = TRUE)
      rads <- Arg(h.mean)
      if(rads<0){rads <- rads + 2*pi}
      averageHue[l] <- rads/(2*pi)
    }
    h.x <- cos(averageHue*2*pi)
    h.y <- sin(averageHue*2*pi)
    h.coords <- cbind(h.x,h.y)
    #print('h coordinates')
    #print(head(h.coords))
    d <- as.matrix(dist(h.coords)) #figure out color in hue space and take difference between them (delta)
    #print('str of d')
    #print(str(d))
    #figure how many genes overlap row-to-row #this is where the bug is, calculating set overlap is slightly wrong
    overlap <- (1-is.na(x))%*%t(1-is.na(x))
    o.diag <- diag(overlap)
    o.square <- sqrt(o.diag%*%t(o.diag))
    overlap <- 1-(overlap/o.square)
    #print("pritning overlap dist")
    #print(overlap)
    #combine these two metrics
    #print('d matrix')
    #print(d[1:5,1:5])
    #print('overlap matrix')
    #print(overlap[1:5,1:5])
    #add distance overlap and color overlap (instead of multiply?????)
    newdist <- weight*d+overlap
    return(as.dist(newdist))
  }
  
  #initiate list of plots
  plots <- list()
  filelist <- list()  
  
  #this loop glues the small matrixes together
  for(i in names(gsea_results[[1]])){
    if (i == "orig_matrix") next
    else {                              #for each experiment V, pull object Vx 
      thisAnalysis <- gsea_results[[i]]
      theBigMatrix <- list() #every time I make a small matrix for each experiment (k), glue them together here
      plots[[i]] <- list()
      #this loop extracts top 30genes from each gmt.file and their top 30 associated pathways, and puts them into a clustered matrix
      for(k in names(gsea_results)){
        print('starting loop with:')
        print(i)
        print(k)
        
        orig_matrix = gsea_results[[k]]$orig_matrix #pull out and separate for appropriate k every time
        theseResults <- gsea_results[[k]][[i]]$Results
        theseResults <- theseResults[theseResults$ES > 0,] #because nmf only returns positives, we dont want negative enrichment scores.
        #write a if, next loop that throws out the garbage (if there is no pathways, don't draw plot)
        if (is.na(theseResults[1,1]) == T) {
          print("These results are EMPTY!!!")
          next
        }
        theseResults <- theseResults[1:min(9, nrow(theseResults)),]
        print("these results")
        print(head(theseResults))#Selects up to top 30 enriched pathways and their stats by normalized enrichment score
        theEdge = unique(unlist(theseResults$leadingEdge)) #Takes all the leading edge genes from theseResults
        theTopGenes = theEdge[order(match(theEdge, orig_matrix$ENTREZID))][1:min(90, length(theEdge))] #Returns the row index of theEdge vs orig_matrix, which is sorted by log2FoldChange (decreasing)
        #print('theTopGenes:')
        #print(length(theTopGenes)) #checks length, should be 30
        agg_leading_edge = orig_matrix[orig_matrix$ENTREZID %in% theTopGenes,c("log2FoldChange","ENTREZID", "SYMBOL")] # ", SYMBOL" Takes the genes corresponding to the index of theTopGenes
        #print('agg_leading_edge:')
        #print(str(agg_leading_edge)) #checks dim of agg_leading_edge, should be 30 x 2
        #print('theseResults$pathway:')
        #print(str(theseResults$pathway)) #checks dim of agg_leading_edge, should be 30 x 2
        #make a small matrix of agg_leading_edge and theseResults$pathways for each k
        if (length(theseResults$pathway) < 1){
          next
        }
        
        theMatrix <- matrix(data = 1,
                            byrow = T,
                            nrow = length(theseResults$pathway),
                            ncol = length(agg_leading_edge$ENTREZID))
        
        print(dim(theMatrix)) #checks dim of theMatrix, should be 30X30 or so
        
        # if (is.null(dim(theMatrix == T))) { #trying to skip plotting empty marices
        #   print("Not enough to plot")
        #   next
        # }
        rownames(theMatrix) <- theseResults$pathway #assign rownames (pathways)
        colnames(theMatrix) <- agg_leading_edge$ENTREZID #assign colnames (genes) #SYMBOL
        
        #for each of the leading edges, clean out theMatrix
        print('theMatrix:')
        print(str(theMatrix))
        
        
        for(j in 1:length(theseResults$leadingEdge)){ 
          #print('j:')
          #print(j)
          tle <- (theseResults$leadingEdge[[j]])
          #tle <- (theseResults$leadingEdge[[j]])##theseResults$leadingEdge ia ENTREZ, theMatrix is now SYMBOL
          theMatrix[j, !(colnames(theMatrix) %in% tle)] <- 0
        }
        #image(t(theMatrix))
        print('theMatrix after j loop')
        print(str(theMatrix))
        
        #rowSD changed to eowMeans etc
        #take standard deviation of the rows and columns
        theMatrix_row_std = rowMeans(theMatrix) # create a vector with all rowSds, keeps the order of rows in theMatrix
        theMatrix_col_std = colMeans(theMatrix)
        
        #print(theMatrix_row_std)
        #print(theMatrix_col_std)
        #remove all rows and cols that have less stdev that 0.01 to avid empty space in plot
        theMatrix <- theMatrix[!(theMatrix_row_std == 0), !(theMatrix_col_std == 0), drop = FALSE] #drop = FALSE makes sure to return a matrix even if not convenient for r
        
        print("str of Matrix before merge")
        print(str(theMatrix))
        
        #somewhox .x and .y is added to genenames in this merge, because genes show up in more than one V, once this is fixed 129 should work
        if(typeof(theBigMatrix)=="list"){
          #first iteration
          theBigMatrix <- theMatrix
          
        }else{
          theBigMatrix <- merge(theBigMatrix,theMatrix,by = "row.names",all=TRUE)
          rownames(theBigMatrix) <- theBigMatrix$Row.names #adding pathway names into row names
          theBigMatrix <- theBigMatrix[,!(names(theBigMatrix) %in% "Row.names")] #removing column with pathways names
          theBigMatrix <- as.matrix(theBigMatrix)
          theBigMatrix[!is.finite(theBigMatrix)] <- 0 #removing NA/NAN/Inf from matrix before plotting
          #image(t(theBigMatrix))
          redundantX <- grep(pattern = ".x", colnames(theBigMatrix), value = T)
          redundantX <- gsub(pattern = "\\.x", replacement = "", x = redundantX)
          for (z in redundantX){
            indX <- grep(colnames(theBigMatrix), pattern = z) #find the col indices of .x and .y
            theBigMatrix[which(theBigMatrix[,indX[2]] == 1), #"if there is a 1 value in the ".y", make it a 1 in the ".x"
                         indX[1]] = 1
            theBigMatrix = theBigMatrix[, -indX[2]] #remove the .y column
            colnames(theBigMatrix)[indX[1]] = z #remove the .x by replacing it with original z
          }
          
        }
      }
      
      print('theBigMatrix at end of loop:')
      print(str(theBigMatrix))
      print('dimensions of theBigMatrix')
      print(dim(theBigMatrix))
      
      
    }
    if (is.null(dim(theBigMatrix == T))) { #trying to skip plotting empty marices
         print("Not enough to plot")
         next
      }
    
  #}
  
    
    print('k is done looping over V1, V2 etc. . . ')
    print('the matrix of genes and genesets has been constructed:')
    print(str(theBigMatrix))
    
    print('theBigMatrix is clean')
    image(theBigMatrix)
   
    #This forloop determines RGB triplet for each gene, builds corresponding matrices and plots a raster
    rbgmatrix <- matrix(data = "#FFFFFF",
                        ncol = length(colnames(theBigMatrix)), 
                        nrow = length(rownames(theBigMatrix))) 
    
    colnames(rbgmatrix) = colnames(theBigMatrix) 
    rownames(rbgmatrix) = rownames(theBigMatrix) 
    
    print("rbg matrix initialized")
    
    for(setID in 1:nrow(theBigMatrix)){ 
      for(geneID in 1:ncol(theBigMatrix)){ 
        tbmo <- theBigMatrix[setID, geneID] 
        if (tbmo == 0){
          rbgmatrix[setID, geneID] <- "#000000"
        }
        else {e <- colnames(theBigMatrix)[geneID] 
        r <- gsea_results$V1$orig_matrix[e,"log2FoldChange"]  #gsea_results$V1$orig_matrix[which(gsea_results$V1$orig_matrix$SYMBOL == e)
                                         #,"log2FoldChange"] 
        g <- gsea_results$V2$orig_matrix[e,"log2FoldChange"]  #gsea_results$V2$orig_matrix[which(gsea_results$V2$orig_matrix$SYMBOL == e)
                                         #,"log2FoldChange"]
        b <- gsea_results$V3$orig_matrix[e,"log2FoldChange"] #gsea_results$V3$orig_matrix[which(gsea_results$V3$orig_matrix$SYMBOL == e)
                                         #,"log2FoldChange"]
        rbgmatrix[setID, geneID] <- rgb(r, g, b)
        }
      }
    }
    print('head of rbgmatrix')
    print(head(rbgmatrix))
    
    ###make shadow hsv matrix that looks like rgbmatrix but instead of #000000 it has hsv.h value, so that we can cluster by color
    hsvmatrix <- matrix(data = 0,
                        ncol = length(colnames(theBigMatrix)), 
                        nrow = length(rownames(theBigMatrix))) 
    
    colnames(hsvmatrix) = colnames(theBigMatrix) 
    rownames(hsvmatrix) = rownames(theBigMatrix) 
    
    print("HSV matrix initialized")
    
    for(setID in 1:nrow(theBigMatrix)){ 
      for(geneID in 1:ncol(theBigMatrix)){ 
        tbmo <- theBigMatrix[setID, geneID] 
        if (tbmo == 0){
          hsvmatrix[setID, geneID] <- NA
        }
        else {e <- colnames(theBigMatrix)[geneID] 
        r <- gsea_results$V1$orig_matrix[e,"log2FoldChange"] # gsea_results$V1$orig_matrix[which(gsea_results$V1$orig_matrix$SYMBOL == e)
                                         #,"log2FoldChange"] 
        g <- gsea_results$V2$orig_matrix[e,"log2FoldChange"]  #gsea_results$V2$orig_matrix[which(gsea_results$V2$orig_matrix$SYMBOL == e)
                                         #,"log2FoldChange"]
        b <- gsea_results$V3$orig_matrix[e,"log2FoldChange"]  #gsea_results$V3$orig_matrix[which(gsea_results$V3$orig_matrix$SYMBOL == e)
                                         #,"log2FoldChange"]
        hsvmatrix[setID, geneID] <- rgb(r, g, b)
        hsv <- t(ColToHsv(hsvmatrix[setID, geneID]))
        hsvmatrix[setID, geneID] <- hsv[,1]
        }
      }
    }
    class(hsvmatrix) <- "numeric"
    print('str of hsvmatrix')
    print(str(hsvmatrix))
    print('clustering hsv matrix now')
    
    ### clustering hsv matrix based on both color and pathway membership
    if (nrow(hsvmatrix) > 2) {
      
      print(hsvmatrix[1:6,1:6])
      heatmap((hsvmatrix),Rowv = NA, Colv = NA,main = "before clustering")
      
      rowcluster <- hclust(HSV.ROW.distfun(hsvmatrix,3))
      colcluster <- hclust(HSV.ROW.distfun(t(hsvmatrix),1.5))
      hsvmatrixOrdered <- hsvmatrix[rowcluster$order, colcluster$order] #ordered rows and columns by rowClust$order
      
      heatmap((hsvmatrix),Rowv = as.dendrogram(rowcluster), Colv =  as.dendrogram(colcluster),main = "after clustering 3, 1.5")
      
      plot(rowcluster)
      plot(colcluster)
    } else { 
      warning("not enough to cluster")
      next
    }
    print('hsvmatrixordered')
    print(str(hsvmatrixOrdered))
    
    #order rbgmatrix rows by hsvmatrixordered rows
    rbgmatrix <- rbgmatrix[rownames(hsvmatrixOrdered),colnames(hsvmatrixOrdered),drop=FALSE]
    
    
    #make tall table of x, y, color for ggplot, use plotly 
    tallmatrix <- as.character(rbgmatrix)
    tallmatrix <- data.frame(data = tallmatrix,
                             x = rep(1:dim(rbgmatrix)[2], each = dim(rbgmatrix)[1]),
                             y = rep(1:dim(rbgmatrix)[1], times = dim(rbgmatrix)[2]),
                             stringsAsFactors = F)
    
    
    xlabels <- c(colnames(rbgmatrix))
    xlabels <- RatENT2RatSYM(xlabels) #add symbols here instead of before
    ylabels <- c(rownames(rbgmatrix))
    
    print(paste0(i, ".pdf"))
    filelist[[i]] <- paste0(i, ".pdf")
    
    pdf(paste0(i, ".pdf"),
        width = ncol(rbgmatrix),
        height= nrow(rbgmatrix))
    
    library(ggplot2)
    
    plots[[i]] <- ggplot(tallmatrix) +
      geom_raster(aes(x = x,
                      y= y,
                      fill = data)) +
      scale_fill_identity() +
      scale_y_continuous(labels = ylabels, breaks = 1:length(ylabels)) +
      scale_x_continuous(labels = xlabels, breaks = 1:length(xlabels)) +
      theme(axis.text.x = element_text(size = 25, angle = 45, vjust = 0.5, hjust=1),
            axis.text.y = element_text(size = 25)) +
      ggtitle(paste0(i))
    
    print(plots[[i]])
    dev.off()
    }
  
  return(filelist)
}






