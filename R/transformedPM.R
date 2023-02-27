#' Transformed Pattern Markers protocol designed to identify the most highly weighted outputs from NMF, and extract only the genes which strongly associate with a particular pattern, or a linear combination of patterns. 
#' @export
#' @import ClusterProfiler
#' @import data.table
#' @param W is a matrix with 3 vectors, the direct un-normalized output of the preceding NMF computation
#' @return nothing
#' @examples 
#' tpm_LLB_NMF <- transformedPM(LLB_NMF)

transformedPM <- function(W){
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