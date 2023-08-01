#' Extract differential expression results from limma model
#'
#' @param fit Limma model fit object
#' @param var_interest Variable of interest 
#' @param var_level Specific level of var_interest to extract 
#' @param numerator Numerator level for comparison
#' @param denominator Denominator level for comparison
#' @param return_tidy Logical to return tidy tibble (default TRUE)
#' @param save_output Logical to save output file  
#' @param output_filename Custom output filename 
#'
#' @return DE results, tidy tibble or full table
#'
#' @examples
#' fit <- rnaseq_to_limma(metadata, counts, design)
#' results <- limma_to_dge(fit, "condition", "treated")
#' 
#' @importFrom limma makeContrasts contrasts.fit topTable
#' @importFrom dplyr select rename
#' @export


limma_to_dge <- function(fit, 
                         var_interest, 
                         var_level = NULL,
                         numerator = NULL,
                         denominator = NULL,
                         return_tidy = TRUE,
                         save_output = NULL,
                         output_filename = NULL) {
  
  #Argument validation
  if(!is.null(var_level) & (!is.null(numerator) | !is.null(denominator))){
    stop("Please specify either 1) var_level OR 2) numerator AND denominator.")
  }
  
  if(sum(is.null(numerator), is.null(denominator)) == 1){
    stop("Please specify numerator AND denominator.")
  }
  
  #Define contrast matrix and coef_name for continuous variables
  if(is.null(var_level) & is.null(numerator) & is.null(denominator)){
    coef_name <- var_interest
    contrast.matrix <- parse(text = paste("makeContrasts(",
                                          coef_name, "=",
                                          var_interest, 
                                          ", levels = colnames(fit$coefficients))") ) %>%
      eval
  }
  
  #Define contrast matrix and coef_name for a coefficient in the model matrix
  if(!is.null(var_level)){
    coef_name <- paste(var_interest, var_level, sep="_")
    contrast.matrix <- parse(text = paste("makeContrasts(",
                                          coef_name, "=",
                                          paste0(var_interest, var_level), 
                                          ", levels = colnames(fit$coefficients))") ) %>%
      eval
  }
  
  #Define contrast matrix and coef_name for a new contrast
  if(!is.null(numerator) & !is.null(denominator)){
    coef_name <- paste(var_interest, numerator, "vs", denominator, sep="_")
    
    if((paste0(var_interest, numerator) %in% colnames(fit$coefficients))&
       (paste0(var_interest, denominator) %in% colnames(fit$coefficients))
    ){
      contrast.matrix <- parse(text = paste("makeContrasts(",
                                            coef_name, "=",
                                            paste0(var_interest, numerator), " - ", 
                                            paste0(var_interest, denominator),
                                            ", levels = colnames(fit$coefficients))") ) %>%
        eval  
    }
    
    
    
    if((paste0(var_interest, numerator) %in% colnames(fit$coefficients))&
       !(paste0(var_interest, denominator) %in% colnames(fit$coefficients))
    ){
      warning("*** Denominator not found in results! *** 
            This could either be beacuse the denominator is the reference level or because the level does not exist.
            Returning the results for the numerator level and assuming the reference level is the demoninator.
            The var_level argument returns the results for a coefficient without generating this warning message.")
      contrast.matrix <- parse(text = paste("makeContrasts(",
                                            coef_name, "=",
                                            paste0(var_interest, numerator),
                                            ", levels = colnames(fit$coefficients))") ) %>%
        eval
    }
    
  }
  
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2,robust = T)
  
  results <- topTable(fit2, number = Inf)
  
  if(return_tidy){
    results <- results %>% 
      dplyr::select(gene_id, logFC,pvalue =  P.Value, fdr = adj.P.Val) 
  }
  
  # Save output
  save <- coalesce(save_output, getOption("save.output"))
  if(save) {
    if(is.null(output_filename)) {
      output_filename <- paste0(coef_name, "limma.dge.Rds")
    }
    saveRDS(results, file.path("output", output_filename))
  }
  
  
  return(results)
  
}