# contains the possible scoring methods for NMF's W

# Description
fc <- function(W){
  #fc <- data.frame(matrix(nrow = nrow(W), ncol = ncol(W)))
  #colnames(fc) <- names(W)
  log2fc <- matrix(nrow = nrow(W), ncol = ncol(W))

  for (sample in 1:ncol(W)) {
    #print(sample)
    curCol <- W[,sample]
    otherCol <- W[,-c(sample)]
    otherCol <- apply(otherCol, 1, mean)
      
    #fc[,sample] <- curCol/otherCol
    log2fc[,sample] <- log2(curCol/otherCol)
  #rownames(fc) <- rownames(W)
  #rownames(log2fc) <- rownames(W)
  }
  #log2fc <- as.matrix(log2fc)
  rownames(log2fc) <- rownames(W)
  
  return(log2fc)
}


# Description
patternMarker <- function(W){
  normedMatrix <- t(apply(W, 1, function(row) row / max(row)))
  
  unitVector <- function(n, length)
  {
    vec <- rep(0, length)
    vec[n] <- 1
    return(vec)
  }
  
  markerScores <- sapply(1:ncol(normedMatrix), function(patternIndex)
    apply(normedMatrix, 1, function(row)
    {
      lp <- unitVector(patternIndex, ncol(normedMatrix))
      return(sqrt(sum((row-lp)^2)))
    })
  )
}


#description
percentScore <- function(W){
  score <- t(apply(W, 1, function(row) row/sum(row)))
  #score <- as.data.frame(score)
  row.names(score) <- row.names(W)
  
  return(score)
}


# description
diffAvg <- function(W){
  scored <- matrix(ncol = ncol(W), nrow = nrow(W))
  for(col in 1:ncol(W)){
    if (ncol(W) > 2){
    scored[,col] <- W[,col] - rowMeans(W[,-col])
    }
    else{
      scored[,col] <- W[,col] - (W[,-col]) # can't rowMeans one col.
    }
  }
  rownames(scored) <- rownames(W)
  return(scored)
}