 
#' Extract shrunken DESeq2 results (apeglm)
#'
#' @param dds DESeqDataSet after differential expression analysis  
#' @param var_interest Variable of interest
#' @param numerator Numerator value for variable
#' @param denominator Denominator value for variable
#' @param use_adaptiveshrink Logical, use adaptive prior (default TRUE)
#'
#' @return Shrunken DESeq2 results as tidy table
#'
#' @examples 
#' dds <- DESeq(dds)
#' results <- deseq2_to_apeglm_dge(dds, var_interest="condition",
#'                                 numerator=1, denominator=0)
#' 
#' @importFrom DESeq2 results lfcShrink 
#' @importFrom tibble rownames_to_column
#' @export

deseq2_to_apeglm_dge <- function(dds,
                                 var_interest,
                                 numerator = NULL,
                                 denominator = NULL,
                                 use_adaptiveshrink = TRUE) {
  
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
  
  # Extract & format results 
  res <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    dplyr::select(gene_id, log2FoldChange, pvalue, padj) %>%
    dplyr::rename(logFC = log2FoldChange,
           fdr = padj)
  
  return(res)
  
}