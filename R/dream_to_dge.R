#' Convert Dream results to tidy table
#' 
#' @param dream_model Fitted Dream model object
#' @param contrast Contrast coefficients
#'
#' @return Tidy tibble with DE results 
#'
#' @examples
#' library(variancePartition)
#' metadata <- tibble(sample = c("s1","s2","s3"),
#'                    type = c("treated","untreated","treated"))
#' counts <- tibble(gene_id = c("ENSG1", "ENSG2"),  
#'                  s1 = c(10, 20),
#'                  s2 = c(15, 50),
#'                  s3 = c(20, 100))
#' design <- model.matrix(~ type, metadata)
#' dream_model <- dream(counts[,2:3], design)
#' results <- dream_to_dge(dream_model, contrast = c(-1,1))
#'
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column 
#' @export

dream_to_dge <- function(dream_model, contrast) {
  
  fit <- contrasts.fit(dream_model, contrast)
  results <- topTable(fit)
  
  results <- results %>%
    rownames_to_column("gene_id") %>%
    select(gene_id, logFC, P.Value, adj.P.Val) %>%
    rename(logFC = logFC, 
           pvalue = P.Value,
           fdr = adj.P.Val)
  
  return(results)
  
}