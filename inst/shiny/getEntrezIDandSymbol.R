
# takes in genelist, which should be the named list of stats from
# analysisres$clusterlist$genelist, and will return a dataframe with columns stats, entrezid
# and symbol
genelistWithEntrezandSymbol <- function(genelist, species){
  # select correct species. Would like more flexibility in this.
  if(species== "Mus musculus"){
    require(org.Mm.eg.db)
    db = org.Mm.eg.db
    species = "Mus musculus"
  }
  else if(species == "Homo sapiens"){
    require(org.Hs.eg.db)
    db = org.Hs.eg.db
  }
  else if(species == "Rattus norvegicus"){
    require(org.Rn.eg.db)
    db = org.Rn.eg.db
  }
  
  
  # get SYMBOL from ENTREZID
  df <- clusterProfiler::bitr(names(genelist),
                                  fromType = "ENTREZID", toType = "SYMBOL",
                                  OrgDb = db)
  # add statistic
  df$stat <- genelist
  
  return(df)
}
