#' Convert Ensembl IDs to HGNC symbols
#'
#' @param results Table with Ensembl gene IDs
#' 
#' @return Table with HGNC symbols
#'
#' @examples 
#' results <- deseq2_to_dge(dds)
#' results <- ensembl_to_hgnc(results)
#'
#' @importFrom dplyr left_join select rename
#' @export

ensembl_to_hgnc <- function(results) {
  
  # Confirm IDs are Ensembl
  id_type <- detect_id_type(results$gene_id)
  if(id_type != "ensembl") {
    stop("IDs are not Ensembl")
  }
  
  # Load mapping
  ensembl_map <- system.file("data", "ensembl_to_symbol.Rds", package = "baton", mustWork = T) %>%
    readRDS()
  
  # Convert IDs
  results <- left_join(results, ensembl_map, by=c("gene_id"= "ensembl_gene_id")) 
  
  return(results)
  
}

# Detect ID type
detect_id_type <- function(ids) {
  if(all(grepl("^ENSG", ids))) {
    return("ensembl")
  } else {
    stop("Not Ensembl IDs")
  }
}