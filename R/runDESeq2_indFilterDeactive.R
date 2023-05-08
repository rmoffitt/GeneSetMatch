# DESeq2
#' This function takes raw count data and performs Differential Expression Analysis of sequence count data based on negative binomial distribution
#' default: independentFiltering = FALSE, cooksCut = .99 quantile of the F(p, m-p)
#' #' We reference DESeq2 for this protocol. Citation:  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq
#' data with DESeq2 Genome Biology 15(12):550 (2014)
#' @export
#' @import biomaRt
#' @import dplyr
#' @import readr
#' @import gplots
#' @import ggplot2
#' @import DESeq2
#' @param count.dat is the raw count data 
#' @param conditions is the sample information for DE
#' @return Results of DE analysis (pval, adjusted pval, logFC, avgExpr, stat)
#' @examples runDESeq2_indFilterDeactive(snyder_readcounts, full_Snyder$sampInfo$`hnf4a status')


runDESeq2_indFilterDeactive <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(DESeq2)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  # Initialize DESeq2 object
  colData <- data.frame(conditions)
  rownames(colData) <- colnames(count.dat)
  d_deseq2 <- DESeqDataSetFromMatrix(countData=count.dat,
                                     colData=colData,
                                     design=~conditions)
  
  # Differential expression analysis
  d_deseq2 <- DESeq(d_deseq2)
  res  <- results(d_deseq2,  independentFiltering = FALSE)
  res$padj <- ifelse(is.na(res$padj), 1, res$padj)
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(res),
                       logFC=res$log2FoldChange,
                       avgExpr=res$baseMean,
                       pval=res$pvalue,
                       padj=res$padj,
                       stat=res$stat)
  
  return(list(counts=d_deseq2, test=res, summary=sum_df, time=clock))
}