#This script installes all of the required R packages to run GeneSetMatch
#make the function check if a package is already installed 
installCran <- function(package_name){
  install.packages(package_name)
}

installBioc <- function(package_name){
  BiocManager::install(package_name)
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
installCran("dplyr")
installCran("FactoMineR")
installCran("factoextra")
installCran("gplots")
installCran("ggplot2")
installCran("ggpubr")
installCran("ggrepel")
installCran("ggraph")
installCran("gridExtra")
installCran("heatmap.plus")
installCran("msigdbr")
installCran("NMF")
installCran("pheatmap")
installCran("pdist")
installCran("readr")
installCran("RColorBrewer")
installCran("scales")
installCran("shiny")
installCran("shinydashboard")
installCran("shinyjqui")
installCran("shinyWidgets")
installCran("shinycssloaders")
installCran("tibble")
installCran("circlize")
installCran("shinyalert")

# #Bioconductor packages
installBioc("biomaRt")
installBioc("ClusterProfiler")
installBioc("DESeq2")
installBioc("enrichplot")
installBioc("vsn")
installBioc("apeglm")
installBioc("org.Hs.eg.db")
installBioc("org.Rn.eg.db")
installBioc("org.Mm.eg.db")
installBioc("ReactomePA")
installBioc("qusage")

# #GitHub packages
devtools::install_github("kevinblighe/EnhancedVolcano")
devtools::install_github("jokergoo/ComplexHeatmap")