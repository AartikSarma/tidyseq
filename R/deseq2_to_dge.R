
#' Convert DESeq2 results to tidy table
#'
#' @param dds DESeq2 object after DE analysis 
#' @param contrast Contrast coefficients
#'
#' @return Tidy tibble with shrunken DE results
#'
#' @examples
#' library(DESeq2)
#' metadata <- tibble(sample = c("s1","s2","s3"),
#'                    type = c("treated","untreated","treated"))
#' counts <- tibble(gene_id = c("ENSG1","ENSG2"),
#'                  s1 = c(10,20), 
#'                  s2 = c(15,50),
#'                  s3 = c(20,100))
#' dds <- DESeqDataSetFromMatrix(counts[,2:3], metadata, ~ type) 
#' dds <- DESeq(dds)
#' results <- deseq2_to_dge(dds, contrast=c(-1,1))  
#'
#' @importFrom DESeq2 results
#' @importFrom tibble rownames_to_column
#' @export

deseq2_to_dge <- function(dds, contrast) {
  
  # Extract results for given contrast
  res <- results(dds, contrast=contrast, alpha=0.05)
  
  # Convert to tidy format
  res <- res %>% 
    rownames_to_column("gene_id") %>% # gene IDs from rownames
    select(gene_id, log2FoldChange, pvalue, padj) %>%
    rename(logFC = log2FoldChange, # rename columns
           fdr = padj)
  
  return(res)
  
}