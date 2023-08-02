#' Differential expression analysis with DESeq2
#'
#' Perform differential expression analysis on count data
#' 
#' @param metadata Metadata table 
#' @param counts Count data
#' @param design Model design formula 
#'
#' @return DESeq2 result object
#'
#' @examples
#' metadata <- tibble(sample = c("s1","s2","s3"), 
#'                    type = c("treated","untreated","treated"))
#' counts <- tibble(gene_id = c("ENSG1","ENSG2"),
#'                  s1 = c(10,20),
#'                  s2 = c(15,50),
#'                  s3 = c(20,100)) 
#' design <- ~ type 
#' de_results <- rnaseq_to_deseq2(metadata, counts, design)
#'
#' @importFrom DESeq2 DESeqDataSet DESeq DESeqDataSetFromMatrix
#' @export

rnaseq_to_deseq2 <- function(metadata, counts, design) {
  
  validate_data(metadata, counts)
  
  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                colData = metadata,
                                design = as.formula(design))
  
  # Differential expression analysis
  dds <- DESeq2::DESeq(dds, parallel = T)
  
  return(dds)
  
}