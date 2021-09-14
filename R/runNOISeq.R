# NOISeq
#' This function performs analysis on RNA-seq expression data or other similar kind of data. Exploratory plots to evualuate saturation, count distribution, expression per chromosome, type of detected features, features length, etc. Differential expression between two experimental conditions with no parametric assumptions.
#' @export
#' @import biomaRt
#' @import dplyr
#' @import readr
#' @import gplots
#' @import ggplot2
#' @import NOISeq
#' @import snow
#' @import samr
#' @param count.dat is the raw count data 
#' @param conditions is the sample information for DE
#' @return Results of DE analysis (pval, adjusted pval, logFC, theta, prob.DE)
#' @examples runNOISeq(snyder_readcounts, full_Snyder$sampInfo$`hnf4a status')
#' 
#default setting
runNOISeq <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(NOISeq)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat  <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                  rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),] 
  ## Start the clock
  ptm <- proc.time()
  
  #Define factors
  myfactors = data.frame(Group = conditions, GroupRun=paste0(conditions, "_", rep(1, length(conditions))), 
                         Run = paste0("R_", rep(1, length(conditions))))
  
  #convert data into NOISeq object
  mydata <- readData(data = count.dat, factors = myfactors)
  
  
  #Test DE with NOISeqBIO and TMM normalization
  mynoiseqbio = noiseqbio(mydata, norm = "tmm", factor = "Group", filter = 0) #filter=0 means no filtering is performed.
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  #Results
  res <-  as.data.frame(mynoiseqbio@results[[1]])
  res$id <- rownames(res)
  res$q <- 1-res$prob  #This is not equivalent to 1 - adjusted p-valeu. We used it to select genes with 
  #higher probability of DE, which is equivalent to the smallest q here
  res$q[is.na(res$q)] <- 1 
  
  
  
  ## Create summary table
  sum_df <- data.frame(ID=res$id,
                       logFC=-1*res$log2FC,
                       avgExpr=(res[,1]+res[,2])/2,
                       prob.DE=res$prob,
                       theta = -1*res$theta,
                       pval = res$q, # this is not the classical p-value
                       padj = res$q)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  return(list(counts=count.dat, test=res, summary=sum_df, time=clock))
}