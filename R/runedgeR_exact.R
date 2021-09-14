# EdgeR Exact
#' This function takes raw count data and performs Differential Expression analysis on specified condition
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
#' @examples runedgeR_exact(snyder_readcounts, full_Snyder$sampInfo$`hnf4a status')

runedgeR_exact <- function(count.dat, conditions, colors='', path='',specification='')
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
  
  ## estimate common dispersion
  y <- estimateCommonDisp(y)
  
  ## estimate gene specific dispersion
  y <- estimateTagwiseDisp(y)
  
  ## DE test
  et <- exactTest(y)                                # Test
  de_table <- topTags(et,n=nrow(et),adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## visualize samples after normalization
  # createMDS(y_norm,colors=colors,path=path,
  #          specification=paste(specification,'edgeR_norm',sep=""))
  
  ##Visualize Dispersion
  plotname <-paste(path,'Dispersion_edgeR_classic_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  plotBCV(y)
  dev.off()
  
  ## P-values and adjusted p-values
  pval <- et[[1]][,3]
  padj <- p.adjust(pval,method="BH") 
  # We checked that padj gives exactly the same result as FDR in topTags
  
  ## Calculate average expression
  avgExpr <- 2^et$table$logCPM * mean(y[[2]]$lib.size)/(10^6)
  
  ## Return rank gene in toptable
  rank <- match(rownames(count.dat),rownames(de_table))
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(et),
                       logFC=et$table$logFC,
                       avgExpr=avgExpr,
                       pval=pval,
                       padj=padj,
                       rank=rank)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  
  ## Automatic classification by method
  clas <- decideTestsDGE(et, p=0.05, adjust="BH")   
  
  return(list(counts=y, test=et, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}