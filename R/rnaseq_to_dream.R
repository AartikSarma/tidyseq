#' Differential expression analysis with Dream
#'
#' Perform differential expression analysis with mixed effects models on RNA-seq count data using Dream
#'
#' @param metadata Metadata table
#' @param counts Count data table
#' @param design Model design formula
#' @param ddf Method to estimate degrees of freedom
#' 
#' @return Limma model fit object
#'
#' @examples
#' metadata <- readRDS("metadata.rds")
#' counts <- readRDS("counts.rds")
#' design <- ~ condition + batch
#' model <- rnaseq_to_dream(metadata, counts, design)
#' 
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @importFrom limma  voom lmFit eBayes  
#' @importFrom variancePartition voomWithDreamWeights dream
#' @importFrom BiocParallel SnowParam
#' @importFrom dplyr filter
#' @export

rnaseq_to_dream <- function(metadata, counts,
                            design,
                            filterByExpr.design = ~1, 
                            ddf = "Kenward-Roger"
                            ) {
  
  validate_data(metadata, counts)
  
  param = SnowParam(6, "SOCK", progressbar=TRUE)
  
  # Create design matrix
  mm <- model.matrix(design, metadata)
  mm.filter <- model.matrix(filterByExpr.design, metadata)
  
  # Create DGEList object
  dge <- edgeR::DGEList(counts, # count data
                 genes = rownames(counts), # gene ids
                 samples = metadata, # sample data
                 remove.zeros = T)
  
  # Normalize 
  dge <- edgeR::calcNormFactors(dge)
  
  # Filter low counts
  keep <- edgeR::filterByExpr(dge, design = mm.filter)
  
  dge <- dge[keep, ]
  #Voom 
  v <- variancePartition::voomWithDreamWeights(counts = dge,
                            formula = design,
                            data = metadata,
                            BPPARAM = param,
                            plot = T)
  #Dream 
  fit <- variancePartition::dream(v,
               formula = design,
               data = metadata, 
               BPPARAM = param,
               ddf = ddf
  ) %>% variancePartition::eBayes(robust = T)
  
  return(fit)
  
}
