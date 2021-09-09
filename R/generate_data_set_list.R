#' generate_data_set_list
#'
#' @export

generate_data_set_list <- function() {
  data_set_list <- data.frame(labels=c(
    "Desired Display Name here"),
    variablenames =c("Actual name of RData object when you load it"))
  
  save(list = c("data_set_list"),
       file = "./data/data_set_list.RData")
  return(NULL)
}
       
