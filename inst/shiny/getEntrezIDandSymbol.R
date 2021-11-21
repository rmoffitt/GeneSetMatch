# takes in OrigMatrix, which should be the matrix from 
# analysisres$clusterlist$OrigMatrix, which contains a stat and ENTREZID column, 
# and will return a matrix with columns named stats, entrezid and symbol
genelistWithEntrezandSymbol <- function(origMatrix, species){
  # select correct species. Would like more flexibility in this.
  if(species == "Mus musculus"){
    require(org.Mm.eg.db)
    db = org.Mm.eg.db
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
  df <- clusterProfiler::bitr(origMatrix$ENTREZID,
                                  fromType = "ENTREZID", toType = "SYMBOL",
                                  OrgDb = db)
  # add statistic
  df$stat <- origMatrix$stat
  
  return(df)
}
