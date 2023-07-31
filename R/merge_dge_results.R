#' Merge multiple DGE results tables
#'
#' @param ... DGE results tables  
#' @param suffixes Optional suffixes
#'
#' @return Merged DGE results
#'
#' @examples
#' table1 <- deseq2_to_dge(dds1)
#' table2 <- deseq2_to_dge(dds2) 
#' merged <- merge_dge_results(table1, table2)
#'
#' @importFrom plyr join_all  
#' @importFrom stringr str_c
#' @export

merge_dge_results <- function(..., suffixes = seq_along(...)) {
  
  # Get number of tables
  n <- length(list(...))
  
  # Handle suffixes
  if(missing(suffixes)) {
    suffixes <- seq_len(n)
  } else {
    if(!is.character(suffixes)) {
      stop("suffixes must be a character vector")
    }
    if(length(suffixes) != n) {
      stop("suffixes length must match number of tables")
    }
  }
  
  # Suffix tables
  results <- list(...)
  for(i in seq_len(n)) {
    suffix <- suffixes[[i]]
    results[[i]] <- results[[i]] %>%
      rename_at(vars(-gene_id), ~str_c(.x, suffix))
  }
  
  # Merge results
  merged <- plyr::join_all(results, type = "full", by = "gene_id")
  
  return(merged)
  
}