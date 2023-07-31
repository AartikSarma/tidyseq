#' Extract differential expression results from limma model
#'
#' @param fit Limma model fit object  
#' @param var_interest Variable of interest
#' @param numerator Numerator level
#' @param denominator Denominator level  
#' @param save_output Logical to save output
#' @param output_filename Custom filename for output
#'
#' @return Tibble with DE results 
#'
#' @examples
#' fit <- rnaseq_to_limma(metadata, counts, design)
#' results <- limma_to_dge(fit, "condition", 1, 0)
#'
#' @importFrom limma contrasts.fit topTable
#' @importFrom dplyr select rename 
#' @export

limma_to_dge <- function(fit, 
                         var_interest,
                         numerator = NULL,
                         denominator = NULL,
                         save_output = NULL,
                         output_filename = NULL) {
  
  # Construct coefficient name
  coef_name <- make_coef_name(var_interest, numerator, denominator) 
  
  # Extract results
  results <- extract_results(fit, coef_name)
  
  # Format results
  results <- format_results(results)
  
  # Save output
  save <- coalesce(save_output, getOption("save.outputs"))
  if(save) {
    if(is.null(output_filename)) {
      output_filename <- paste0(coef_name, ".Rds")
    }
    saveRDS(results, file.path("output", output_filename))
  }
  
  return(results)
  
}

# Internal functions
make_coef_name <- function(var, num, den) {
  if(is.null(num) & is.null(den)) {
    return(var)
  } else {  
    return(paste(var, num, "vs", den, sep="_"))
  }
}

extract_results <- function(fit, coef) {
  contrast.matrix <- makeContrasts(contrasts=coef,
                                   levels=colnames(fit$coefficients))
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2, robust=T)
  
  return(topTable(fit2, number=Inf))
}

format_results <- function(res) {
  res %>%
    select(gene_id, logFC, P.Value, adj.P.Val) %>% 
    rename(logFC = logFC,
           pvalue = P.Value, 
           fdr = adj.P.Val)
}