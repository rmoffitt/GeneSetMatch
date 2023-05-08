#' Updating formatting of experiment files 
#'
#' \code{Moff2Summarized} will convert a standard list format experiment file into a SummarizedExperiment
#'
#' @param object An object of class 'list' with 4 items: $ex - the expression matrix or dataframe, $sampInfo - the colData, $featInfo - the rowData, and $metadata
#' 
#' @return Returns an object of class 'SummarizedExperiment' that contains all info from inital experiment contained in those 4 list headings,
#' see ?\code{SummarizedExperiment} for help with manipulating new object
#'

#' @export
#' @import SummarizedExperiment

Moff2Summarized <- function(object) {
   if(class(object) == "SummarizedExperiment"){
      warning("Object is already a Summarized Experiment, nothing is done")
   }
   else if(class(object) != "list"){
      warning(paste0("Object is of class '", class(object), "' not a list, please check your input and formatting,
              Required formatting is:
                        object = list(ex = expression matrix or df,
                                      sampInfo = df of column info
                                      featInfo = df of row info,
                                      metadata = list of metadata)"))
   }
   else if(!exists("ex",where = object)){
      warning("Expression data does not exist, nothing to convert")
   }
   else{
      new_obj <- SummarizedExperiment(assays = object$ex,
                                      rowData = object$featInfo,
                                      colData = object$sampInfo,
                                      metadata = object$metadata)
      return(new_obj)
   }
}
