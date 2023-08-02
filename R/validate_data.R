#' Validate metadata and expression data
#'
#' Check that metadata and expression formats are compatible
#'
#' @param metadata Metadata table
#' @param expression Gene expression data
#'
#' @return NULL if validation passes, stops with error otherwise  
#'
#' @examples
#' metadata <- readRDS("metadata.rds")
#' expression <- readRDS("expression.rds")
#' validate_data(metadata, expression)
#'
#' @export

validate_data <- function(metadata, expression) {
  
  # Check expression data is numeric
  if(!all(sapply(expression, is.numeric))) {
    stop("Expression data contains non-numeric columns")
  }
  
  # Check expression has rownames
  if(is.null(rownames(expression))) {
    stop("Expression data must have rownames that are gene identifiers")
  }
  
  # Check metadata colnames match expression colnames
  meta_samples <- rownames(metadata)
  exp_samples <- colnames(expression)
  if(!all(exp_samples %in% meta_samples)) {
    stop("Expression column names not in metadata row names") 
  }
  
  # Check metadata rownames match expression colnames
  if(!all(meta_samples %in% exp_samples)) {
    stop("Metadata row names not in expression column names")
  }
  
  # Check equal
  if(!identical(meta_samples, exp_samples)) {
    stop("Metadata and expression samples not identical")
  }
  
  # If reach here, valid
  return(NULL)
  
}