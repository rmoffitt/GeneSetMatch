#For Debugging Purposes
NEWpreparePathwaysAndStats <- function(pathways, stats, minSize, maxSize, gseaParam, scoreType){
  # Error if pathways is not a list
  if (!is.list(pathways)) {
    stop("pathways should be a list with each element containing names of the stats argument")
  }
  
  # Error if stats is not named
  if (is.null(names(stats))) {
    stop("stats should be named")
  }
  
  # Warning message for ties in stats
  ties <- sum(duplicated(stats[stats != 0]))
  if (ties != 0) {
    warning("There are ties in the preranked stats (",
            paste(round(ties * 100 / length(stats), digits = 2)),
            "% of the list).\n",
            "The order of those tied genes will be arbitrary, which may produce unexpected results.")
  }
  
  # Warning message for duplicate gene names
  if (any(duplicated(names(stats)))) {
    warning("There are duplicate gene names, fgsea may produce unexpected results.")
  }
  
  if (all(stats > 0) & scoreType == "std"){
    warning("All values in the stats vector are greater than zero and scoreType is \"std\", ",
            "maybe you should switch to scoreType = \"pos\".")
  }
  
  stats <- sort(stats, decreasing=TRUE)
  stats <- abs(stats) ^ gseaParam
  
  
  minSize <- max(minSize, 1)
  
  pathwaysFiltered <- lapply(pathways, function(p) { unique(na.omit(match(p, names(stats)))) })
  pathwaysSizes <- sapply(pathwaysFiltered, length)
  
  toKeep <- which(minSize <= pathwaysSizes & pathwaysSizes <= maxSize)
  
  pathwaysFiltered <- pathwaysFiltered[toKeep]
  pathwaysSizes <- pathwaysSizes[toKeep]
  
  list(filtered=pathwaysFiltered,
       sizes=pathwaysSizes,
       stats=stats)
}