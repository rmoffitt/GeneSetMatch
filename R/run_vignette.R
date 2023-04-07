#' Wrapper function to facilitate calling running of vignettes when package is installed, not cloned.
#'
#' @export
#'
#' @param format A character string of length one that specifies which vignette to build. Must be one of "bulk", "NMF", or "scRNA".
#'
#' @import rmarkdown


run_vignette <- function(format = c("bulk","NMF","scRNA")){
   if(format == "bulk"){ ## Est run time 5 min
      rmarkdown::render(system.file("analysis/RunMe_bulk.Rmd", package = "GeneSetMatch"),
                        output_dir = getwd())
   }
   if(format == "NMF"){ ## Est run time 15 min
      rmarkdown::render(system.file("analysis/RunMe_NMF.Rmd", package = "GeneSetMatch"),
                        output_dir = getwd())
   }
   if(format == "scRNA"){ ## Est run time 15 min
      rmarkdown::render(system.file("analysis/RunMe_SingleCell.Rmd", package = "GeneSetMatch"),
                        output_dir = getwd())
   }
   return(NULL)
}
