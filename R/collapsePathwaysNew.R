#For Debugging Purposes
collapsePathwaysNew <- function(fgRes,
                             pathways,
                             stats,
                             pval.threshold=0.05,
                             nperm=10/pval.threshold,
                             gseaParam=1) {
  universe <- names(stats)
  
  pathways <- pathways[fgRes$pathway]
  pathways <- lapply(pathways, intersect, universe)
  
  print("fgRes")
  print(fgRes)
  
  parentPathways <- setNames(rep(NA, length(pathways)), names(pathways))
 
  for (i in seq_along(pathways)) {
    p <- names(pathways)[i]
    if (!is.na(parentPathways[p])) {
      next
    }

    pathwaysToCheck <- setdiff(names(which(is.na(parentPathways))), p) 
    print("pathways to check")
    print(pathwaysToCheck)
    
    pathwaysUp <- fgRes[pathways %in% pathwaysToCheck & fgRes$NES >= 0, "pathway"]#fgRes[which(fgRes$pathway %in% pathways & fgRes$ES >= 0),"pathway"] #returns only enriched pathways as opposed to "[pathway %fin% pathwaysToCheck & ES >= 0][, pathway]"
    pathwaysDown <- fgRes[pathways %in% pathwaysToCheck & fgRes$NES < 0, "pathway"]#fgRes[which(fgRes$pathway %in% pathways & fgRes$ES < 0),"pathway"] ## maybe I'm wrong, look intro replacing p with pathways or something/anything with length >1, causes problems downstream in plotting because it trims everything to just 1
    print(paste0("PathwaysUp: ", head(pathwaysUp))) #added print command
    print(paste0("PathwaysDwn: ", head(pathwaysDown)))
    
    if (length(pathwaysToCheck) == 0) {
      break
    }
    
    minPval <- setNames(rep(1, length(pathwaysToCheck)), pathwaysToCheck)
    
    u1 <- setdiff(universe, pathways[[p]])
    
    fgResUp1 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u1],
                               nperm=nperm, maxSize=length(u1)-1, nproc=1,
                               gseaParam=gseaParam, scoreType = "pos")
    fgResDown1 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u1],
                                 nperm=nperm, maxSize=length(u1)-1, nproc=1,
                                 gseaParam=gseaParam, scoreType = "neg")
    fgRes1 <- rbindlist(list(fgResUp1, fgResDown1), use.names = TRUE)
    
    minPval[fgRes1$pathway] <- pmin(minPval[fgRes1$pathway], fgRes1$pval)
    
    u2 <- pathways[[p]]
    
    fgResUp2 <- fgseaSimple(pathways = pathways[pathwaysUp], stats=stats[u2],
                               nperm=nperm, maxSize=length(u2)-1, nproc=1,
                               gseaParam=gseaParam, scoreType = "pos")
    fgResDown2 <- fgseaSimple(pathways = pathways[pathwaysDown], stats=stats[u2],
                                 nperm=nperm, maxSize=length(u2)-1, nproc=1,
                                 gseaParam=gseaParam, scoreType = "neg")
    fgRes2 <- rbindlist(list(fgResUp2, fgResDown2), use.names = TRUE)
    
    minPval[fgRes2$pathway] <- pmin(minPval[fgRes2$pathway], fgRes2$pval)
    
    parentPathways[names(which(minPval > pval.threshold))] <- p
  }
  
  return(list(mainPathways=names(which(is.na(parentPathways))),
              parentPathways=parentPathways))
}