#' Validate metadata and gene expression data
#'
#' Confirm metadata and gene expression data formats are compatible for analysis
#'
#' @param metadata Tibble with sample metadata 
#' @param gene_expression Tibble with gene expression data
#' @param sample_name Name of sample ID column in metadata (default 'sampleid')
#' @param gene_id Name of gene ID column in expression (default 'ensembl_gene_id')
#' 
#' @return TRUE if validation passed, FALSE otherwise
#'
#' @examples
#' metadata <- tibble(sampleid = c("s1", "s2", "s3"))
#' expression <- tibble(ensembl_gene_id = c("ENSG1", "ENSG2"),
#'                       s1 = c(1,2),
#'                       s2 = c(2,1), 
#'                       s3 = c(3,5))
#' validate_data(metadata, expression) # Returns TRUE
#' 
#' @export
validate_data <- function(metadata, gene_expression, 
                          sample_name = "sampleid",
                          gene_id = "ensembl_gene_id") {
  
  # Check columns
  meta_samples <- unique(metadata[[sample_name]])
  exp_samples <- colnames(gene_expression)[-1]
  
  if(!all(meta_samples %in% exp_samples)) {
    return(FALSE)
  }
  
  if(nrow(metadata) != length(exp_samples)-1) {
    return(FALSE) 
  }
  
  if(!is.character(gene_expression[[1]])) {
    return(FALSE)
  } 
  
  if(sum(sapply(gene_expression[,2:ncol(gene_expression)], is.numeric)) != ncol(gene_expression)-1) {
    return(FALSE)
  }
  
  return(TRUE)
  
}