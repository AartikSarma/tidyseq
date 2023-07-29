#' Differential expression analysis with limma and voom
#'
#' Perform differential expression analysis on count data using limma + voom
#'
#' @param metadata Metadata table  
#' @param counts Count data
#' @param design Model design formula 
#'
#' @return Fitted limma model object
#'
#' @examples
#' metadata <- tibble(sample = c("s1","s2","s3"),
#'                    type = c("treated","untreated","treated"))
#' counts <- tibble(gene_id = c("ENSG1","ENSG2"),
#'                  s1 = c(10,20), 
#'                  s2 = c(15,50),
#'                  s3 = c(20,100))
#' design <- ~ type
#' limma_model <- rnaseq_to_limma(metadata, counts, design)  
#'
#' @importFrom limma voom  lmFit eBayes
#' @importFrom edgeR DGEList
#' @export

rnaseq_to_limma <- function(metadata, counts, design) {
  
  # Voom transformation
  v <- voom(counts[,2:ncol(counts)], design, plot=FALSE)
  
  # Create DGEList
  dge <- DGEList(v$E, genes=counts$gene_id)
  
  # Linear model fit
  design <- model.matrix(design, data=metadata) 
  dge <- lmFit(dge, design)
  
  # Empirical Bayes statistics
  dge <- eBayes(dge)
  
  return(dge)
  
}