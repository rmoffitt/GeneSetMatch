#' Write a gmt file.
#' @export
#' @import dplyr
#' @import data.table
#' @param list_of_sets A named list where each element in the list is an output of msigdbr
#' @return nothing
#' @examples
#' C1_df = msigdbr(species = "Mus musculus", category = "C1", subcategory = "")
#' C2_df = msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
#' list_of_sets = list(C1 = C1_df,
#'                    C2 = C2_df)
#' convert_msigdbr_obj_to_gmt_file(list_of_sets)              

convert_msigdbr_obj_to_gmt_file <- function(list_of_sets) {
  file_name_list <- c()
  for(i in names(list_of_sets)){
    to_parse = list_of_sets[[i]]
    m_t2g  = to_parse %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
    my_aggregated_list <- aggregate(m_t2g$entrez_gene, by = list(gs = m_t2g$gs_name), FUN = function(x){
      paste(x, sep = "\t", collapse = "\t")})
    my_aggregated_list <- data.frame(gs = my_aggregated_list$gs, 
                                     url = "lucie.com", 
                                     genes = my_aggregated_list$x)
    my_var <- paste0("./","my", i,".gmt")
    write.table(my_aggregated_list, 
                file = my_var, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    print(paste0("wrote .gmt for ",i))
    file_name_list[i] <- my_var
  }
  return(file_name_list)
  }

