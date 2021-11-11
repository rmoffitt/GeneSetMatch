# prepares genesets to be compatible with ODIS.GSEA


# issue -- update to allow species selection and category/subcategory selection, as well as user submitted pathways
prepareGenesets <- function(species, selectPathList){
  library(msigdbr)
  list_of_sets = list()
  for (i in selectPathList) {
    print(i)
    list_of_sets[[i]] = assign(paste(i, "_df", sep = ""), msigdbr(species = species, category = i))
  }
  
  #C1_df = msigdbr(species = species, category = "C1")
  #C2_df = msigdbr(species = species, category = "C2")
  #C2KEGG_df = msigdbr(species = species, category = "C2")
  #C3_df = msigdbr(species = species, category = "C3")
  #C4_df = msigdbr(species = species, category = "C4")
  #C5_df = msigdbr(species = species, category = "C5")
  #C6_df = msigdbr(species = species, category = "C6")
  #C7_df = msigdbr(species = species, category = "C7")
  #C8_df = msigdbr(species = species, category = "C8")

  #list_of_sets = list(C1 = C1_df,
                      #C2 = C2_df,
                      #C2_KEGG = C2KEGG_df,
                      #C3 = C3_df,
                      #C4 = C4_df,
                      #C5 = C5_df,
                      #C6 = C6_df,
                      #C7 = C7_df,
                      #C8 = C8_df
                      #)

  source("~/ODIS2/R/convert_msigdbr_obj_to_gmt_file.R") # issue
  list_of_gmt_files <- convert_msigdbr_obj_to_gmt_file(list_of_sets)
  
  return(list_of_gmt_files)
}


# be sure to match species above
ODISGSEA_helper <- function(gene_list, gmtList, pval){ # issue? need to think about how this will work
  GSEA_outputs <- list()
    for(i in gmtList){
    print(i)
    source("~/ODIS2/R/ODISGSEA.R") # issue
    res <- ODISGSEA(gene_list = gene_list, theGoFile = i, pval = pval)
    GSEA_outputs[["gene_list"]][[i]] <- res ## need to structure outputs correctly
  }
  
  return(GSEA_outputs)
}
