#' Updating formatting of experiment files
#'
#' \code{Moff2Seurat} will convert a standard list format experiment file into a Seurat Object
#'
#' @param object An object of class 'list' with at least 3 items: $ex - the expression matrix or dataframe, $sampInfo - the colData, $featInfo - the rowData, and $metadata (optional)
#' 
#' @param project A character string of length one to be given as the project name to the Seurat object, defaults to the name of object
#' 
#' @param new_names Either a numeric index or character string of length one of which column in $featInfo to use as rownames
#'  
#' @return Returns an object of class 'SeuratObj' that contains all info from inital experiment contained in those 4 list headings,
#' see ?\code{Seurat} for help with manipulating new object
#'

#' @export
#' @import Seurat

Moff2Seurat <- function(object, project = NULL, new_names = 1) {
   if(class(object) == "Seurat"){
      warning("Object is already a Seurat Object, nothing is done")
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
      rownames(object$ex) = object$featInfo[,new_names]
      if(is.null(project)){
         project = as.character(deparse(substitute(object)))
      }
      new_obj <- CreateSeuratObject(counts = object$ex,
                                    meta.data = object$sampInfo,
                                    project = project)
      if(exists("metadata",where = object)){ new_obj@misc = object$metadata }
      return(new_obj)
   }
}
