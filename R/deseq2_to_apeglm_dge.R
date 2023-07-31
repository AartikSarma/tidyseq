#' Extract shrunken DESeq2 results (apeglm)
#'
#' @param dds DESeqDataSet after DE analysis  
#' @param var_interest Variable of interest
#' @param numerator Numerator value (default NULL)
#' @param denominator Denominator value (default NULL)
#' @param use_adaptiveshrink Logical for adaptive prior (default TRUE) 
#' @param save_output Logical to save results (default NULL) 
#' @param output_filename Custom filename for results (default NULL)
#'
#' @return Shrunken DESeq2 results as tidy table
#'
#' @examples
#' \dontrun{
#' options(save.output = TRUE)
#' results <- deseq2_to_apeglm_dge(dds, save_output = FALSE)  
#' }
#'
#' @importFrom DESeq2 results lfcShrink
#' @importFrom tibble rownames_to_column
#' @export

deseq2_to_apeglm_dge <- function(dds, 
                                 var_interest, 
                                 numerator = NULL,
                                 denominator = NULL,
                                 use_adaptiveshrink = TRUE,
                                 save_output = NULL,
                                 output_filename = NULL) {
  
  # Construct coefficient name
  ifelse((is.null(numerator) & is.null(denominator)), 
         coef_name <- var_interest, 
         coef_name <- paste(var_interest, numerator, "vs", denominator, sep="_")
  )
  
  
  # Apply apeglm shrinkage
  res <- lfcShrink(dds, 
                   coef=coef_name, 
                   type="apeglm",
                   apeAdapt =use_adaptiveshrink)
  
  # Format results
  res <- format_results(res)
  
  # Determine if saving outputs
  save <- coalesce(save_output, getOption("save.output"))

  if(save) {
    
    # Set output filename if not provided
    if(is.null(output_filename)) {
      output_filename <- paste0(coef_name, ".apeglm.dge.Rds") 
    }
    
    # Save results
    saveRDS(res, file.path("output", output_filename))
  }
  
  return(res)
  
}

# Helper function  
format_results <- function(res) {
  res %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    select(gene_id, log2FoldChange, pvalue, padj) %>% 
    rename(logFC = log2FoldChange,  
           fdr = padj)
}