# QUASISeq QL
#' This function identifies differentially expressed genes in RNA-seq count data using quasi-Poisson or quasi-negative binomial models.
#' We reference QuasiSeq in this protocol. Citation: Lund S, Qu L, Bhasavanich K, Marschner I, R Core Team (2022). _QuasiSeq: Analyzing RNA
#' Sequencing Count Tables Using Quasi-Likelihood_. R package version 1.0-11-0,
#' <https://CRAN.R-project.org/package=QuasiSeq>.
#' @export
#' @import biomaRt
#' @import dplyr
#' @import readr
#' @import gplots
#' @import ggplot2
#' @import QuasiSeq
#' @import edgeR
#' @import snow
#' @import samr
#' @param count.dat is the raw count data 
#' @param conditions is the sample information for DE
#' @return Results of DE analysis (pval, adjusted pval, logFC, stat)
#' @examples runQuasiSeq_QL(snyder_readcounts, full_Snyder$sampInfo$`hnf4a status')
#' 


runQuasiSeq_QL <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(QuasiSeq)
  require(edgeR)
  
  conditions <- as.factor(conditions)
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                             rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  ## Start the clock
  ptm <- proc.time()
  
  ##Estimating dispersion using edgeR package (trended dispersion estimate from edgeR GLM)
  y <- DGEList(counts=count, group=conditions)
  y <- calcNormFactors(y)
  design <- model.matrix(~conditions)
  est.nb.disp <- estimateGLMTrendedDisp(y, design)$trended.dispersion
  
  ### Create list of designs describing model under null and alternative hypotheses
  design.list <- vector("list",2)
  design.list[[1]] <- model.matrix(~as.factor(conditions)) #This also could have just been ``trt''.
  design.list[[2]] <- rep(1,length(conditions))
  log.offset <- log(apply(count,2,quantile,.75))
  
  ### Analyze using QL, QLShrink and QLSpline methods applied to quasi-Poisson model
  ### applied to quasi-negative binomial model
  fit <- QL.fit(count, design.list, log.offset=log.offset, nb.disp=est.nb.disp, Model="NegBin", 
                print.progress=FALSE, NBdisp="trend")
  results <- QL.results(fit, Dispersion="Deviance", Plot=TRUE)
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID   =rownames(count),
                       logFC=fit$coefficients[,2], 
                       pval = as.vector(results$P.values$QL),
                       padj = as.vector(results$Q.values$QL),
                       stat = as.vector(results$F.stat$QL))
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  
  return(list(counts=y, test=fit,  summary=sum_df, time=clock))
}