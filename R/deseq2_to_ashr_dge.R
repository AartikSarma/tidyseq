 
#' Extract shrunken DESeq2 results (ashr)
#'
#' @param dds DESeqDataSet after differential expression analysis
#' @param var_interest Variable of interest
#' @param numerator Numerator value for variable  
#' @param denominator Denominator value for variable
#'
#' @return Shrunken DESeq2 results as tidy table
#'
#' @examples
#' dds <- DESeq(dds)  
#' results <- deseq2_to_ashr_dge(dds, var_interest="condition", 
#'                               numerator=1, denominator=0)
#'
#' @importFrom DESeq2 results lfcShrink  
#' @importFrom tibble rownames_to_column
#' @export

deseq2_to_ashr_dge <- function(dds, 
                               var_interest, 
                               numerator,
                               denominator) {
  
  # Construct contrast
  contrast <- setNames(object = c(numerator, denominator),
                       nm = levels(dds[[var_interest]]))
  
  # Apply ashr shrinkage 
  res <- lfcShrink(dds, coef=contrast, type="ashr")
  
  # Extract & format results
  res <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    select(gene_id, log2FoldChange, pvalue, padj) %>%
    rename(logFC = log2FoldChange,  
           fdr = padj)
  
  return(res)
  
}