#This script installes all of the required R packages to run GeneSetMatch
#make the function check if a package is already installed 
installCran <- function(package_name){
  if(!require(package_name,character.only = TRUE)){
      install.packages(package_name)
    }
}

installBioc <- function(package_name){
  if(!require(package_name,character.only = TRUE)){
    BiocManager::install(package_name)
  }
}

#may need to run this
options(repos = "https://cloud.r-project.org")
options(install.packages.compile.from.source="interactive")

#cran packages
installCran("BiocManager")
installCran("BBmisc")
installCran("colourpicker")
installCran("data.table")
installCran("devtools") # to install from github
installCran("DT")
installCran("DescTools")
installCran("dplyr")
installCran("FactoMineR")
installCran("factoextra")
installCran("gplots")
installCran("ggplot2")
installCran("ggpubr")
installCran("ggrepel")
installCran("ggraph")
installCran("gridExtra")
installCran("msigdbr")
installCran("matrixStats")
installCran("NMF")
installCran("pheatmap")
installCran("pdist")
installCran("plotly")
installCran("readr")
installCran("raster")
installCran("RColorBrewer")
installCran("scales")
installCran("Seurat")
installCran("shiny")
installCran("shinydashboard")
installCran("shinyjqui")
installCran("shinyWidgets")
installCran("shinycssloaders")
installCran("snow")
installCran("samr")
installCran("sp")
installCran("stringr")
installCran("tibble")
installCran("circlize")
installCran("shinyalert")
installCran("QuasiSeq")
installCran("xfun")

# #Bioconductor packages
installBioc("annotate")
installBioc("baySeq")
installBioc("biomaRt")
installBioc("BiocParallel")
installBioc("clusterProfiler")
installBioc("DESeq2")
installBioc("enrichplot")
installBioc("edgeR")
installBioc("fgsea")
installBioc("gage")
installBioc("genefilter")
installBioc("impute")
installBioc("limma")
installBioc("NOISeq")
installBioc("ROTS")
installBioc("vsn")
installBioc("apeglm")
installBioc("org.Hs.eg.db")
installBioc("org.Rn.eg.db")
installBioc("org.Mm.eg.db")
installBioc("rat2302.db")
installBioc("ReactomePA")
installBioc("SummarizedExperiment")
installBioc("qusage")

# #GitHub packages
devtools::install_github("kevinblighe/EnhancedVolcano")
devtools::install_github("jokergoo/ComplexHeatmap")

