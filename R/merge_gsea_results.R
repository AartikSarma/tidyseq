#' Merge multiple GSEA result tables
#'
#' Bind rows of multiple GSEA result tables into one table
#' 
#' @param ... List of GSEA result data frames
#' @param gsea_result_names Custom names for each GSEA result  
#'
#' @return Merged data frame of all GSEA results
#'
#' @examples
#' results1 <- dge_to_gsea(results1)
#' results2 <- dge_to_gsea(results2)
#' merged <- merge_gsea_results(results1, results2, 
#'                              gsea_result_names = c("Treatment1", "Treatment2"))
#'
#' @importFrom dplyr bind_rows mutate
#' @export

merge_gsea_results <- function(..., gsea_result_names = NULL) {
  
  # Validate custom names if provided
  if(!is.null(gsea_result_names)) {
    if(length(gsea_result_names) != length(list(...))) {
      stop("gsea_result_names must match number of results")
    }
  }
  
  # Default names if none provided
  if(is.null(gsea_result_names)) {
    gsea_result_names <- paste0("Result", seq_along(list(...)))
  }
  
  # Bind rows 
  merged <- bind_rows(list(...), .id = "gsea_result") %>%
    mutate(gsea_result = as.numeric(gsea_result)) %>% # convert id to numeric
    mutate(gsea_result_name = gsea_result_names[gsea_result]) # add names
  
  return(merged)
  
}