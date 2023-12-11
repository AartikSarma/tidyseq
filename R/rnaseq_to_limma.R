#' Differential expression analysis with limma-voom  
#'
#' Perform differential expression analysis on RNA-seq count data using limma-voom
#'
#' @param metadata Metadata table
#' @param counts Count data table
#' @param design Model design formula
#' 
#' @return Limma model fit object
#'
#' @examples
#' metadata <- readRDS("metadata.rds")
#' counts <- readRDS("counts.rds")
#' design <- ~ condition + batch
#' model <- rnaseq_to_limma(metadata, counts, design)
#' 
#' @importFrom edgeR DGEList filterByExpr calcNormFactors voomLmFit
#' @importFrom limma  voom lmFit eBayes  
#' @importFrom dplyr filter
#' @export

rnaseq_to_limma <- function(metadata, counts, design,
                            filterByExpr.design = ~1) {
  
  validate_data(metadata, counts)
  
  # Create design matrix
  mm <- model.matrix(design, metadata)
  mm.filter <- model.matrix(filterByExpr.design, metadata)
  
  # Create DGEList object
  dge <- DGEList(counts, # count data
                 genes = rownames(counts), # gene ids
                 samples = metadata, # sample data
                 remove.zeros = T)
  
  # Normalize 
  dge <- calcNormFactors(dge)
  
  # Filter low counts
  keep <- filterByExpr(dge, design = mm.filter)
  
  #voomLmFit
  fit <- voomLmFit(
    counts = dge[keep,], 
    design = mm, 
    plot = T
  ) %>% eBayes(robust = T)
   

  return(fit)
  
}