#' Convert Dream results to tidy table
#' 
#' @param dream_model Fitted Dream model object
#' @param contrast Contrast coefficients
#'
#' @return Tidy tibble with DE results 
#'
#' @importFrom dplyr rename 
#' @export

dream_to_dge <- function(dream_model, coef) {
  
  results <- variancePartition::topTable(dream_model ,coef = coef, number = Inf) %>% 
    dplyr::select(gene_id = genes, 
                  logFC = logFC, 
                  pvalue = P.Value,
                  fdr = adj.P.Val)
  
  return(results)
  
}