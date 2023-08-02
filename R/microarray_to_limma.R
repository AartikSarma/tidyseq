#' Differential expression analysis for microarray data
#'
#' Perform DE analysis on microarray data using limma
#'
#' @param metadata Metadata table
#' @param expression Microarray expression matrix 
#' @param design Model design formula 
#'
#' @return Fitted limma model object
#'
#' @examples
#' metadata <- tibble(sample = c("s1","s2","s3"),
#'                    type = c("treated","untreated","treated"))
#' expression <- matrix(rnorm(9*200), nrow = 9, ncol = 200) 
#' design <- ~ type
#' limma_model <- microarray_to_limma(metadata, expression, design)
#' 
#' @importFrom limma lmFit eBayes
#' @export

microarray_to_limma <- function(metadata, expression, design) {
  validate_data(metadata, expression)
  
  mm <- model.matrix(design, data = metadata)
  fit <- lmFit(expression, mm)
  fit <- eBayes(fit, robust = TRUE)
  
  return(fit)
  
}