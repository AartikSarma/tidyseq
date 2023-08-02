#' Merge multiple result tables into long format
#'
#' Bind rows of multiple result tables into one long table
#'
#' @param ... List of result data frames/tibbles
#' @param result_names Custom names for each result object
#'
#' @return Merged long-format data frame  
#' 
#' @examples
#' results1 <- deseq2_to_dge(dds1)
#' results2 <- deseq2_to_dge(dds2)
#' merged <- merge_results_long(results1, results2,
#'                              result_names = c("Treatment1", "Treatment2"))
#'
#' @importFrom dplyr bind_rows mutate  
#' @export

merge_results_long <- function(..., result_names = NULL) {
  
  # Validate consistent column names
  col_names <- lapply(list(...), colnames)
  if(length(unique(col_names)) != 1) {
    stop("Column names not consistent across objects")
  }
  
  # Validate custom result names if provided
  if(!is.null(result_names)) {
    if(length(result_names) != length(list(...))) {
      stop("result_names must match number of objects")
    }
  }
  
  # Default result names if none provided
  if(is.null(result_names)) {
    result_names <- paste0("Result", seq_along(list(...))) 
  }
  
  # Bind rows into single table
  merged <- bind_rows(list(...), .id = "result") %>%
    mutate(result = as.numeric(result)) %>% # convert id to numeric
    mutate(result_name = result_names[result]) # add names
  
  return(merged)
  
}