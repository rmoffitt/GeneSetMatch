#' Menu of all the GeneSetMatch tools for Differential Expression Analysis 
################################# edgeR  ################################
#' @export
#' @import biomaRt
#' @import dplyr
#' @import readr
#' @import gplots
#' @import ggplot2
#' @import edgeR
#' @import QuasiSeq
#' @import limma
#' @import ROTS
#' @import NOISeq
#' @import baySeq
#' @import snow
#' @import samr
#' @param list_of_sets A named list where each element in the list is an output of msigdbr
#' @return nothing
#' @examples
#' C1_df = msigdbr(species = "Mus musculus", category = "C1", subcategory = "")
#' C2_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
#' list_of_sets = list(C1 = C1_df,
#'                    C2 = C2_df)
#' convert_msigdbr_obj_to_gmt_file(list_of_sets)       



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

##  edgeR GLM  
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

## edgeR robust   
### default (prior.df=10)
runedgeR_rob_pDF10 <- function(count.dat, conditions, colors='', path='',specification='')
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
  
  ## Set up design matrix and run GLM
  design <- model.matrix(~conditions)
  y <- estimateGLMRobustDisp(y, design)   # as in pipeline.pdf document
  
  ## DE test
  dFit <- glmFit(y, design)                         
  dLrt <- glmLRT(dFit,coef=2)                       # Test
  de_table <- topTags(dLrt,n=nrow(dLrt),adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  
  ## Visualize Dispersion
  plotname <-paste(path,'Mean_dispersion_edgeR_rob_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  par(mar=c(5.1,6.1,4.1,2.1))
  plotBCV(y, xlim=c(-5, 12), ylim=c(0, 5), cex.axis=1, cex.lab=1, cex=1)
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
  
  return(list(counts=y, test=dLrt, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

### prior.df=5
runedgeR_rob_pDF05 <- function(count.dat, conditions, colors='', path='',specification='')
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
  
  ## Set up design matrix and run GLM
  design <- model.matrix(~conditions)
  y <- estimateGLMRobustDisp(y, design, prior.df = 5)   # as in pipeline.pdf document
  
  ## DE test
  dFit <- glmFit(y, design)                         
  dLrt <- glmLRT(dFit,coef=2)                       # Test
  de_table <- topTags(dLrt,n=nrow(dLrt),adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  
  ## Visualize Dispersion
  plotname <-paste(path,'Mean_dispersion_edgeR_rob_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  par(mar=c(5.1,6.1,4.1,2.1))
  plotBCV(y, xlim=c(-5,12), ylim=c(0, 5), cex.axis=2, cex.lab=2, cex=1)
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
  
  return(list(counts=y, test=dLrt, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

### (prior.df=20)
runedgeR_rob_pDF20 <- function(count.dat, conditions, colors='', path='',specification='')
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
  
  ## Set up design matrix and run GLM
  design <- model.matrix(~conditions)
  y <- estimateGLMRobustDisp(y, design, prior.df = 20)   # as in pipeline.pdf document
  
  ## DE test
  dFit <- glmFit(y, design)                         
  dLrt <- glmLRT(dFit,coef=2)                       # Test
  de_table <- topTags(dLrt,n=nrow(dLrt),adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  
  ## Visualize Dispersion
  plotname <-paste(path,'Mean_dispersion_edgeR_rob_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  par(mar=c(5.1,6.1,4.1,2.1))
  plotBCV(y, xlim=c(-5, 10), ylim=c(0,5), cex.axis=2, cex.lab=2, cex=1)
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
  
  return(list(counts=y, test=dLrt, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}


### prior.df calculated based on 
runedgeR_rob_auto.pDF <- function(count.dat, conditions, colors='', path='',specification='')
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
  pDF <- getPriorN(y)  #Calculating the appropriate prior degrees of freedom
  
  ## Normalize
  y <- calcNormFactors(y)
  
  ## Set up design matrix and run GLM
  design <- model.matrix(~conditions)
  y <- estimateGLMRobustDisp(y, design, prior.df = pDF)   # as in pipeline.pdf document
  
  ## DE test
  dFit <- glmFit(y, design)                         
  dLrt <- glmLRT(dFit,coef=2)                       # Test
  de_table <- topTags(dLrt,n=nrow(dLrt),adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  
  ## Visualize Dispersion
  plotname <-paste(path,'Mean_dispersion_edgeR_rob_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  par(mar=c(5.1,6.1,4.1,2.1))
  plotBCV(y, xlim=c(-5, 12), ylim=c(0,5), cex.axis=2, cex.lab=2, cex=1)
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
  
  return(list(counts=y, test=dLrt, tt=de_table, de_aut=clas, summary=sum_df, time=clock, calculated.prior.df=pDF))
}


## edgeR quasi-likelihood
runedgeR_ql<-function(count.dat, conditions, colors='', path='',specification='') 
  {
  require(edgeR)
  
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  y <- DGEList(counts=count.dat, group=conditions)
  y <- calcNormFactors(y)
  design <- model.matrix(~conditions)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  dFit <- glmQLFTest(fit)    
  de_table <- topTags(dFit, n=nrow(dFit), adjust.method='BH',
                      sort.by='PValue')             # Test results ordered by p-value
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## P-values and adjusted p-values
  pval <- dFit$table$PValue
  padj <- p.adjust(pval,method="BH")
  # It was checked that padj gives exactly the same result as FDR
  
  ## Calculate average expression
  avgExpr <- 2^dFit$table$logCPM * mean(y[[2]]$lib.size)/(10^6)
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(dFit),
                       logFC=dFit$table$logFC,
                       avgExpr=avgExpr,
                       pval=pval,
                       padj=padj,
                       stat=dFit$table$F)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by method
  clas <- decideTestsDGE(dFit, p=0.05, adjust="BH") 
  
  return(list(counts=y, test=dFit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}


################################# QuasiSeq #################################
## QL
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
## QLShrink
runQuasiSeq_QLShrink <- function(count.dat, conditions, colors='', path='',specification='')
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
                       pval = as.vector(results$P.values$QLShrink),
                       padj = as.vector(results$Q.values$QLShrink),
                       stat = as.vector(results$F.stat$QLShrink))
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  
  return(list(counts=y, test=fit,  summary=sum_df, time=clock))
}
## QLSpline
runQuasiSeq_QLSpline <- function(count.dat, conditions, colors='', path='',specification='')
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
                       pval = as.vector(results$P.values$QLSpline),
                       padj = as.vector(results$Q.values$QLSpline),
                       stat = as.vector(results$F.stat$QLSpline))
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  
  return(list(counts=y, test=fit,  summary=sum_df, time=clock))
}


##################################### DESeq ####################################

runDESeq <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(DESeq)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  ## Initialize new DESeq object
  cds <- newCountDataSet(count.dat, conditions)
  
  ## estimate size factor
  cds <- estimateSizeFactors(cds)
  cds_norm <- cds
  
  ## estimate dispersion parameters
  cds <- estimateDispersions(cds, sharingMode="maximum",
                             method= "pooled", fitType='parametric') #fitType is "parametric" by default 
  
  ## differential expression
  res <- nbinomTest(cds, levels(conditions)[1], levels(conditions)[2])
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  # Visualize Dispersion
  plotname <-paste(path,'Dispersion_DESeq_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  DESeq::plotDispEsts(cds)
  dev.off()
  
  # Histogram unadjusted p-values
  plotname <-paste(path,'Hist_pval_DESeq_',specification,'.png',sep="")
  png(plotname, width=10, height=7, units="in", res=300)
  hist(res$pval,breaks=100,main="")
  dev.off()
  
  
  ## Create summary table
  sum_df <- data.frame(ID=res$id,
                       logFC=res$log2FoldChange,
                       avgExpr=res$baseMean,
                       pval=res$pval,
                       padj=res$padj)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  return(list(counts=cds, test=res, summary=sum_df, time=clock))
}


##################################### DESeq2 ###################################
##independentFiltering = FALSE, cooksCut = .99 quantile of the F(p, m-p)
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

## independentFiltering = FALSE, cooksCut = OFF
runDESeq2_indFilterDeactiveCooksOff <- function(count.dat, conditions, colors='', path='',specification='')
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
  res  <- results(d_deseq2,  independentFiltering = FALSE, cooksCutoff = FALSE)
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

##default: independentFiltering = TRUE, cooksCut = .99 quantile of the F(p, m-p)
runDESeq2_indFilterActive <- function(count.dat, conditions, colors='', path='',specification='')
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
  res  <- results(d_deseq2,  independentFiltering = TRUE)
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

runDESeq2_indFilterActive_2 <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(DESeq2)
  #' This function is similar to runDESeq2_indFilterActive() except that it doesn't 
  #' replace NA p-adj values by 1. It is simply used to compare number of filtered genes
  #' 
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat  <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
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
  res  <- results(d_deseq2,  independentFiltering = TRUE)
  #res$padj <- ifelse(is.na(res$padj), 1, res$padj)
  
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



################################ baySeq ###################################
# Function to return M and A values for baySeq output 
MA.values <- function(cD, samplesA, samplesB, normaliseData = TRUE)
{
  ## Modified version of plotMA.CD taken from baySeq (v.1.12.0) source code
  ## returns M and A values
  if(is.character(samplesA)) {
    Asamps <-  which(as.character(cD@replicates) %in% samplesA)
    if(!all(samplesA %in% cD@replicates))
      Asamps <- c(Asamps, which(colnames(cD@data) %in% samplesA[!(samplesA %in% as.character(cD@replicates))]))
    if(!all(samplesA %in% c(colnames(cD@data), as.character(cD@replicates)))) warning("Some members of 'samplesA' were not found!")
    samplesA <- Asamps
  }
  if(length(samplesA) == 0)
    stop("Can't find any data for sample set A!")
  
  if(is.character(samplesB)) {
    Bsamps <-  which(as.character(cD@replicates) %in% samplesB)
    if(!all(samplesB %in% cD@replicates))
      Bsamps <- c(Bsamps, which(colnames(cD@data) %in% samplesB[!(samplesB %in% as.character(cD@replicates))]))
    if(!all(samplesB %in% c(colnames(cD@data), as.character(cD@replicates)))) warning("Some members of 'samplesB' were not found!")
    samplesB <- Bsamps
  }
  
  if(length(samplesB) == 0)
    stop("Can't find any data for sample set B!")  
  
  if(!inherits(cD, what = "countData"))
    stop("variable 'cD' must be of or descend from class 'countData'")
  
  Adata <- cD@data[,samplesA]
  Bdata <- cD@data[,samplesB]
  
  if(normaliseData) {
    Adata <- Adata /matrix(rep(libsizes(cD)[samplesA],nrow(Adata)),byrow=TRUE,
                           ncol=ncol(Adata)) * mean(libsizes(cD)[c(samplesA, samplesB)])
    Bdata <- Bdata /matrix(rep(libsizes(cD)[samplesB],nrow(Bdata)),byrow=TRUE,
                           ncol=ncol(Bdata)) * mean(libsizes(cD)[c(samplesA, samplesB)])
  }  
  
  if(nrow(seglens(cD)) > 0)
    if(ncol(seglens(cD)) == 1) {
      Adata <- Adata / seglens(cD)[,1]
      Bdata <- Bdata / seglens(cD)[,1]
    } else {
      Adata <- Adata / seglens(cD)[,samplesA]
      Bdata <- Bdata / seglens(cD)[,samplesB]
    }
  
  Adata <- colSums(t(Adata)) / length(samplesA)
  Bdata <- colSums(t(Bdata)) / length(samplesB)
  
  Azeros <- which(Adata == 0)
  Bzeros <- which(Bdata == 0)
  nonzeros <- which(Adata != 0 & Bdata != 0)
  infRatio <- ceiling(max(abs((log2(Adata) - log2(Bdata))[nonzeros]), na.rm = TRUE))
  
  M <- log2(Adata) - log2(Bdata)
  M[Azeros] <- -infRatio - 2
  M[Bzeros] <- infRatio + 2
  
  A <- (log2(Adata) + log2(Bdata)) / 2
  A[Azeros] <- log2(Bdata[Azeros])
  A[Bzeros] <- log2(Adata[Bzeros])
  
  ma <- cbind(-M, A)  # -M to make it consistent with foldchanges of the other methods
  colnames(ma) <- c("M", "A")
  return(ma)
}

runBaySeq<- function(count.dat, conditions, colors='', path='',specification='')
{
  require(baySeq)
  require(snow)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  ## Create cluster
  if(require("parallel")) cl<-makeCluster(8) else cl<-NULL
  
  ## Create count data object
  cd <- new("countData", data=count.dat, replicates=conditions,
            groups=list(NDE=rep(1,length(conditions)),
                        DE=c(rep(1,length(grep(levels(conditions)[1], conditions))),
                       rep(2,length(grep(levels(conditions)[2], conditions))))))
  
  ## Estimate require scaling factors 
  libsizes(cd) <- getLibsizes(cd, estimationType ='edgeR')  #default is quantile
  
  ## Generate priors
  
  cd <- getPriors.NB(cd, cl=cl)
  
  ## Generate likelihoods
  cd = getLikelihoods.NB(cd, cl=cl)
  
  ## Returns posterior likelihoods
  res <- topCounts(cd, group='DE', normaliseData=TRUE, number = nrow(count.dat))
  
  ## Close cluster
  if(!is.null(cl)) stopCluster(cl)
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Calculate MA values
  ma <- MA.values(cd, samplesA=levels(conditions)[1], samplesB=levels(conditions)[2], normaliseData = TRUE)
  
  ## Extract p-values, adjusted p-values and test statistic
  id <- match(row.names(count.dat), row.names(res))
  pval <- res[id,]$FDR.DE    # pval and padj are identical
  padj <- res[id,]$FDR.DE
  stat <- res[id,]$Likelihood
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(count.dat),
                       logFC=ma[,'M'],
                       avgExpr=2^ma[,'A'],
                       pval=pval,
                       padj=padj,
                       stat=stat)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  ## Return results
  return(list(test=cd, tt=res, summary=sum_df, time=clock))
}




################### limma  #######################
## limma + Quantile Normalization
runLimmaQN <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  # Create model matrix and set up dge object
  design <- model.matrix(~conditions)
  
  # Quantile normalization  
  counts.log.dat=log2(count.dat+1)
  counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')
  dat=counts.log.norm.dat
  
  # Test for differential expression
  fit <- lmFit(dat,design)
  fit <- eBayes(fit)
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dat),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr,
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=dat, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

## limma+Voom (least square estimation) 
runLimmaVoom <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  # Create model matrix and set up dge object
  design <- model.matrix(~conditions)
  dgel <- DGEList(counts=count.dat)
  
  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)
  
  # Voom transformation - estimate dispersion
  v <- voom(dgel,design,plot=FALSE)
  
  # Test for differential expression
  fit <- lmFit(v,design, method="ls")  #default
  fit <- eBayes(fit, trend = FALSE, robust = FALSE)  #default
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr * mean(v[[1]]$lib.size) /(10^6),
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=v, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

## limma+Voom (robust estimation) 
runLimmaVoom_Robust <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  # Create model matrix and set up dge object
  design <- model.matrix(~conditions)
  dgel <- DGEList(counts=count.dat)
  
  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)
  
  # Voom transformation - estimate dispersion
  v <- voom(dgel,design,plot=FALSE)
  
  # Test for differential expression
  fit <- lmFit(v,design, method="ls")  #default
  fit <- eBayes(fit, robust = TRUE, trend = FALSE)
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr * mean(v[[1]]$lib.size) /(10^6),
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=v, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

## limma+trended (least square estimation) 
runLimmaTrended <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  # Create model matrix and set up dge object
  design <- model.matrix(~conditions)
  dgel <- DGEList(counts=count.dat)
  
  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)
  
  # CPM transformation - estimate dispersion
  y <- new("EList")
  y$E <- edgeR::cpm(dgel, log = TRUE, prior.count = 3)
  
  # Test for differential expression
  fit <- lmFit(y, design, method="ls")  #default
  fit <- eBayes(fit, trend = TRUE, robust = FALSE)
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr * mean(dgel$samples$lib.size) /(10^6),
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=y$E, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

## limma+trended (robust estimation) 
runLimmaTrended_Robust <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  # Create model matrix and set up dge object
  design <- model.matrix(~conditions)
  dgel <- DGEList(counts=count.dat)
  
  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)
  
  # CPM transformation - estimate dispersion
  y <- new("EList")
  y$E <- edgeR::cpm(dgel, log = TRUE, prior.count = 3)
  
  # Test for differential expression
  fit <- lmFit(y, design, method="ls")  #default
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr * mean(dgel$samples$lib.size) /(10^6),
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=y$E, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

## limma+Voom (WITH quality weights) 
runLimmaVoom_QW <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  ## Start the clock
  ptm <- proc.time()
  
  # Create model matrix and set up dge object
  design <- model.matrix(~conditions)
  dgel <- DGEList(counts=count.dat)
  
  # Calculate normalization factors - TMM
  dgel <- calcNormFactors(dgel)
  
  # Voom transformation with sample quality weights - estimate dispersion
  v <-  voomWithQualityWeights(dgel,design,normalization="none",plot=FALSE)
  
  # Test for differential expression
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr * mean(v[[1]]$lib.size) /(10^6),
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=v, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}

## limma+Vst  
runLimmaVst <- function(count.dat, conditions, colors='', path='',specification='')
{
  require(limma)
  require(DESeq)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  ## Initialize design matrix and new DESeq object
  design <- model.matrix(~conditions)
  cds <- newCountDataSet(count.dat, conditions)
  
  ## Estimate size factors
  cds <- estimateSizeFactors(cds)
  
  ## Estimate dispersion parameters
  cds <- estimateDispersions(cds,method='blind',fitType='local')
  # Options set as in Soneson research
  
  ## Apply variance stabilizing transformation
  # own function that adapts getVarianceStabilizedData function
  vst <- DESeq::getVarianceStabilizedData(cds)
  
  # Limma test for differential expression
  fit <- lmFit(vst,design)
  fit <- eBayes(fit)
  
  # Return top table
  de_table <- topTable(fit,coef=ncol(design),n=nrow(vst),sort.by="none")
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Create summary table
  sum_df <- data.frame(ID=rownames(de_table),
                       logFC=de_table$logFC,
                       avgExpr=2^de_table$AveExpr,
                       pval=de_table$P.Value,
                       padj=de_table$adj.P.Val,
                       stat=de_table$t)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  # Automatic classification by the method
  clas <- decideTests(fit, p=0.05, adjust='BH')
  
  return(list(counts=cds, test=fit, tt=de_table, de_aut=clas, summary=sum_df, time=clock))
}


################################ SAMSeq ###################################
# Calc mean Expr PoissonSeq Output 
meanEx_SAM <- function(count.dat,conditions)
{
  #Note that here Average Expression is calculated based on the non-normalized counts
  #   (vs. meanEx_Poisson which is done on the average counts)
  conditions <- as.factor(conditions)
  Adata <- count.dat[,conditions==levels(conditions)[1]]
  Bdata <- count.dat[,conditions==levels(conditions)[2]]
  
  A_avg <- colSums(t(Adata)) / sum(conditions==levels(conditions)[1])
  B_avg <- colSums(t(Bdata)) / sum(conditions==levels(conditions)[2])
  
  Azeros <- which(A_avg == 0)
  Bzeros <- which(B_avg == 0)
  
  A <- (log2(A_avg ) + log2(B_avg)) / 2
  A[Azeros] <- log2(B_avg[Azeros])
  A[Bzeros] <- log2(A_avg[Bzeros])
  
  return(A)
}
runSAMSeq <- function(count.dat, conditions, colors='', path='',specification='')
{
  require('samr')
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                 rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  ## Set seed
  set.seed(123)
  
  ## SAM Seq
  x <- count.dat
  y <- factor(conditions)
  capture.output({samfit <- SAMseq(x, y, resp.type = "Two class unpaired",  fdr.output=1,
                                   genenames = rownames(count.dat))})
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  ## Average expression
  A <- 2^meanEx_SAM(count.dat, conditions) 
  
  ## FDR
  f.table = as.data.frame(rbind(samfit$siggenes.table$genes.up, samfit$siggenes.table$genes.lo))
  
  # fdr = rep(NA, nrow(count.dat)) #contains NA value
  # fdr[as.numeric(as.character(f.table[, "Gene Name"]))] <- as.numeric(as.character(f.table[, "q-value(%)"]))/100  
  # 
  ## Create summary table
  temp0 <- data.frame(ID = rownames(samfit$samr.obj$x),
                      avgExp = A[rownames(samfit$samr.obj$x)])
  
  temp1 <- data.frame(ID = as.character(f.table$`Gene ID`), 
                      logFC = log2(as.numeric(as.character(f.table$`Fold Change`))),
                      stat = as.numeric(as.character(f.table$`Score(d)`)),
                      pval = as.numeric(as.character(f.table$`q-value(%)`))/100,
                      padj = as.numeric(as.character(f.table$`q-value(%)`))/100)
    
    
  sum_df <- merge(temp0, temp1, by="ID", all=TRUE) 
  sum_df$pval <- ifelse(is.na(sum_df$pval), 1, sum_df$pval)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  
  ## Extract automatic classification
  # temp3 <- as.data.frame(samfit$siggenes.table$genes.up)
  # temp3$clas <- 1
  # temp4 <- as.data.frame(samfit$siggenes.table$genes.lo)
  # temp4$clas <- -1
  # temp5 <- rbind(temp3,temp4)
  # temp5$ID <- rownames(count.dat)[as.numeric(as.character(temp5[,'Gene Name']))]
  # temp5[,'Score(d)'] <- as.numeric(as.character(temp5[,'Score(d)']))
  # temp5[,'q-value(%)'] <- as.numeric(as.character(temp5[,'q-value(%)']))/100
  # temp5[,'Fold Change'] <- as.numeric(as.character(temp5[,'Fold Change']))
  # clas <- temp5[,c('ID','clas','Score(d)','q-value(%)','Fold Change')]
  # 
  return(list(test=samfit, summary=sum_df, time=clock))
}


################################ NOISeq ###################################
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


################################ ROTS ###################################
#default setting .... based on CPM data
runROTS_cpm <- function(count.dat, conditions, colors='', path='',specification='')
  {
  require(ROTS)
  require(limma)
  require(edgeR)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat  <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                  rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),]
  
  ## Start the clock
  ptm <- proc.time()
  
  #Normalization with TMM
  dge <- DGEList(counts = count.dat)
  dge <- edgeR::calcNormFactors(dge)
  
  #calculate counts per millions from edgeR
  cpms <- cpm(dge)
  
  #Test DE based on log expression
  rots <- ROTS(data = cpms, groups = as.numeric(conditions), B = 1000, K = 1000, log = FALSE)
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  #Results
  res <- data.frame(LFC=rots$logfc, p=rots$pvalue, q=rots$FDR)
  res$q[is.na(res$q)] <- 1
  res$id <- rownames(res) 
  
  ## Calculate average expression
  avgCPMExpr <- rowMeans(rots$data)
  avgExpr    <-  avgCPMExpr * mean(dge$samples$lib.size)/(10^6)
  
  ## Create summary table
  sum_df <- data.frame(ID=res$id,
                       logFC=-1*res$LFC,
                       avgExpr=avgExpr,
                       pval = res$p, # this is not the classical p-value
                       padj = res$q)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  return(list(counts=count.dat, test=res, summary=sum_df, time=clock))
}

#ROTS + voom
runROTS_voom <- function(count.dat, conditions, colors='', path='',specification='')
  {
  require(ROTS)
  require(limma)
  require(edgeR)
  
  #Filtering out genes with 0 count across all samples of each condition 
  conditions <- as.factor(conditions)
  count.dat  <- count.dat[which(rowSums(count.dat[,which(conditions == levels(conditions)[1])])>=1 & 
                                  rowSums(count.dat[,which(conditions == levels(conditions)[2])])>=1),] 
  
  ## Start the clock
  ptm <- proc.time()
  
  #Normalization with TMM
  dge <- DGEList(counts = count.dat)
  dge <- edgeR::calcNormFactors(dge)
  
  #Voom transformation
  design <- model.matrix(~conditions)
  vm  <- voom(dge, design = design)
  
  #Test DE based on log expression
  rots <- ROTS(data = vm$E, groups = as.numeric(conditions), B = 1000, K = 1000, log = TRUE)
  
  ## Stop clock
  clock <- proc.time() - ptm
  
  
  #Results
  res <- data.frame(LFC=rots$logfc, p=rots$pvalue, q=rots$FDR)
  res$q[is.na(res$q)] <- 1
  res$id <- rownames(res) 
  
  ## Calculate average expression
  avgCPMExpr <- 2^de_table$AveExpr * mean(v[[1]]$lib.size) /(10^6)
  avgExpr    <-  avgCPMExpr * mean(dge$samples$lib.size)/(10^6)
  
  ## Create summary table
  sum_df <- data.frame(ID=res$id,
                       logFC=-1*res$LFC,
                       avgExpr=avgExpr,
                       pval = res$p, # this is not the classical p-value
                       padj = res$q)
  sum_df$padj <- ifelse(is.na(sum_df$padj), 1, sum_df$padj)
  return(list(counts=count.dat, test=res, summary=sum_df, time=clock))
}



