# PoissonSEQ
#' This function takes raw count data to implement a method for normalization, testing, and false discovery rate estimation for RNA-sequencing data
#' @export
#' @import biomaRt
#' @import dplyr
#' @import readr
#' @import gplots
#' @import ggplot2
#' @import PoissonSeq
#' @import snow
#' @param count.dat is the raw count data 
#' @param conditions is the sample information for DE
#' @return Results of DE analysis (pval, adjusted pval, logFC, avgExpr, stat)
#' @examples runPoissonSeq(snyder_readcounts, full_Snyder$sampInfo$`hnf4a status')
#' 

runPoissonSeq <- function(count.dat, conditions, colors='', path='',specification='',cutoff=0)
{
  require('PoissonSeq')
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  ## Reformat conditions vector for PoissonSeq analysis
  y <- seq.int(length(conditions))
  y[conditions==levels(conditions)[1]] <- 1
  y[conditions==levels(conditions)[2]] <- 2
  y <- as.numeric(y)
  
  ## Create list object for PoissonSeq analysis
  dat <- list(n=count.dat,
              y=y,
              type='twoclass',
              pair=FALSE,
              gname=rownames(count.dat))
  
  ## Set parameters for Poisson Analysis
  para <- list(ct.sum=1, ct.mean=0.01, npermu=500)
  # Permutations set to 500 as in Rapaport,
  # Exclude genes from DE expression analysis according to threshold set in
  #  function call - 0.5 is standard value for average count
  
  ## Call Poisson DE analysis
  res <- PS.Main(dat=dat,para=para)
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Calculate normalized counts
  lib.factors <- PS.Est.Depth(count.dat)
  norm.count.dat <- as.matrix(count.dat) %*% diag(1/lib.factors)
  A <- meanEx_Poisson(norm.count.dat,conditions)
  
  ## Create summary table
  # Data frame with ID and average expression (includes all genes)
  temp1 <- data.frame(ID=rownames(norm.count.dat), avgExpr=2^A)
  # Data frame with other metrics (only includes genes that past the filter)
  temp2 <- data.frame(ID=rownames(res),
                      logFC=res$log.fc,
                      pval=res$pval,
                      padj=res$fdr,        # Checked: fdr is not exactly the same as padj with BH,
                      # but it is close)
                      stat=res$tt)
  sum_df <- merge(temp1,temp2,by='ID',all.x=TRUE)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  return(list(tt=res, summary=sum_df, time=clock))
}