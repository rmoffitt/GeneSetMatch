# EdgeR GLM
#' This function takes raw count data and performs Differential Analysis of sequence count data via a General Linear Modeling analysis on specified condition
#' We reference EdgeR in this protocol. Citation: Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential
#' expression analysis of digital gene expression data. Bioinformatics 26, 139-140
#' @export
#' @import biomaRt
#' @import dplyr
#' @import readr
#' @import gplots
#' @import ggplot2
#' @import edgeR
#' @param count.dat is the raw count data 
#' @param conditions is the sample information for DE
#' @return Results of DE analysis (pval, adjusted pval, logFC, avgExpr, rank)
#' @examples runedgeR_glm(snyder_readcounts, full_Snyder$sampInfo$`hnf4a status')
#' 

runedgeR_glm <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(edgeR)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  ## Initiate DGElist
  y <- DGEList(counts=count.dat, group=conditions)
  
  ## Normalize
  y <- calcNormFactors(y)
  y_norm <- y
  
  ## Set up design matrix and run GLM
  design <- model.matrix(~conditions)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  # p17 manual tagwise dispersion is strongly recommended in multi-factor 
  #   experiment cases
  
  ## DE test
  dFit <- glmFit(y, design)                         
  dLrt <- glmLRT(dFit,coef=2)                       # Test
  de_table <- topTags(dLrt,n=nrow(dLrt),adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  
  ## visualize samples after normalization
  #createMDS(y_norm, colors=colors, path=path,
  #           specification=paste(specification,'edgeR_norm',sep=""))
  # 
  # # Visualize Dispersion
  plotname <-paste(path,'Mean_dispersion_edgeR_glm_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  par(mar=c(5.1,6.1,4.1,2.1))
  plotBCV(y, xlim=c(-5,12), ylim=c(0,5), cex.axis=2, cex.lab=2, cex=1)
  dev.off()
  
  ## P-values and adjusted p-values
  pval <- dLrt$table$PValue
  padj <- p.adjust(pval,method="BH")
  # It was checked that padj gives exactly the same result as FDR
  
  ## Calculate average expression
  avgExpr <- 2^dLrt$table$logCPM * mean(y[[2]]$lib.size)/(10^6)
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(dLrt),
                       logFC=dLrt$table$logFC,
                       avgExpr=avgExpr,
                       pval=pval,
                       padj=padj,
                       stat=dLrt$table$LR)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by method
  clas <- decideTestsDGE(dLrt, p=0.05, adjust="BH") 
  
  return(list(counts=y, test=dLrt, tt=de_table, de_aut=clas, summary=sum_df, time=clock, y_norm=y_norm))
}