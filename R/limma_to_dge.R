
#' Convert limma results to tidy table
#'
#' @param model Fitted limma model object
#' @param contrast Contrast to extract results for
#'
#' @return Tidy tibble with DE results
#'
#' @examples 
#' metadata <- tibble(sample = c("s1","s2","s3"),
#'                    type = c("treated","untreated","treated"))  
#' expression <- matrix(rnorm(9*200), nrow = 9, ncol = 200)
#' design <- model.matrix(~ type, data = metadata) 
#' limma_model <- lmFit(expression, design)
#' limma_model <- eBayes(limma_model)
#'  
#' # Get results for treated vs untreated
#' de_table <- limma_to_dge(limma_model, contrast = c(-1,1,0)) 
#'
#' @importFrom limma contrasts.fit
#' @importFrom tibble rownames_to_column
#' @export

limma_to_dge <- function(model, contrast) {
  
  fit <- contrasts.fit(model, contrast)
  
  results <- topTable(fit)
  
  results <- results %>% 
    rownames_to_column("gene_id") %>%
    select(gene_id, logFC, P.Value, adj.P.Val) %>%
    rename(logFC = logFC,
           pvalue = P.Value,
           fdr = adj.P.Val)
  
  return(results)
  
}