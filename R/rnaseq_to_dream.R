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
#' @importFrom edgeR DGEList filterByExpr calcNormFactors voomLmFit
#' @importFrom limma  voom lmFit eBayes  
#' @importFrom variancePartition voomWithDreamWeights dream
#' @importFrom dplyr filter
#' @export

rnaseq_to_dream <- function(metadata, counts, design,
                            ddf = "Kenward-Roger",
                            filter_low_counts = TRUE) {
  
  validate_data(metadata, counts)
  
  param = SnowParam(6, "SOCK", progressbar=TRUE)
  
  # Create design matrix
  mm <- model.matrix(design, metadata)
  
  # Create DGEList object
  dge <- DGEList(counts, # count data
                 genes = rownames(counts), # gene ids
                 samples = metadata) # sample data
  
  # Normalize 
  dge <- calcNormFactors(dge)
  
  # Filter low counts
  keep <- filterByExpr(dge, design = mm)
  
  v <- voomWithDreamWeights(counts = dge,
                            formula = design,
                            data = metadata,
                            BPPARAM = param)
  #voomLmFit
  fit <- dream(v,
               formula = design,
               data = metadata, 
               BPPARAM = param,
               ddf = ddf
  ) %>% eBayes(robust = T)
  
  return(fit)
  
}

??dream
