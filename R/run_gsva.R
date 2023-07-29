#' Perform gene set variance analysis
#'
#' Run GSVA on count data to assess pathway activity
#'
#' @param counts Count data as tibble
#' @param metadata Sample metadata as tibble
#' @param design Model design formula
#' @param gsets Gene set database as list of gene vectors
#'
#' @return Tibble with sample scores for each gene set 
#'
#' @examples
#' # Load example data
#' counts <- get_example_counts() 
#' metadata <- get_example_metadata()
#' design <- ~ condition
#' gsets <- get_example_genesets()
#' 
#' # Run GSVA
#' gsva_results <- run_gsva(counts, metadata, design, gsets)
#'
#' @importFrom GSVA gsva
#' @importFrom dplyr left_join
#' @export
run_gsva <- function(counts, metadata, design, gsets) {
  
  # Format count data for GSVA
  count_mat <- counts %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  # Transform counts
  count_mat <- vst(count_mat)
  
  # Run GSVA analysis
  gsva_scores <- gsva(count_mat, gsets, 
                      method="gsva", 
                      kcdf="Gaussian", 
                      abs.ranking=FALSE)
  
  # Convert GSVA scores to tidy format
  gsva_tidy <- as_tibble(t(gsva_scores)) %>% 
    rownames_to_column("sample") %>%
    gather(key="geneset", value="score", -sample) %>%
    left_join(metadata, by="sample")
  
  return(gsva_tidy)
  
}