#' Pattern Markers
#' @export
#' @import ClusterProfiler
#' @import data.table
#' @param W is a matrix with 3 vectors
#' @return nothing
#' @examples 
#' PatternMarker <- LLB_NMF

PatternMarker <- function(W){
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