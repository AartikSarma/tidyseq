#' Merge multiple result tables into wide format
#'
#' Bind rows of multiple result tables and pivot to wide format
#' 
#' @param ... List of result data frames/tibbles
#'
#' @return Merged wide-format data frame
#'
#' @examples 
#' results1 <- deseq2_to_dge(dds1)
#' results2 <- deseq2_to_dge(dds2)
#' merged <- merge_results_wide(results1, results2)
#'
#' @importFrom tidyr pivot_wider
#' @export

merge_results_wide <- function(...) {
  
  # Merge into long format 
  merged_long <- merge_results_long(...)
  
  # Pivot to wide format
  merged_wide <- merged_long %>% 
    pivot_wider(
      id_cols = matches("pathway|gene_id"), # columns to keep unpivoted
      names_from = "result", # column to pivot
      values_from = !(matches("gene_id|pathway|result")), # pivot metric columns
      names_sep = ""
    )
  
  return(merged_wide)
  
}