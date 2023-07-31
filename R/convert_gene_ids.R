
#' Convert gene IDs to symbols
#'
#' @param results Tidy DE results table 
#' @param species Species name (eg. "human", "mouse")
#' @param id_type Species name (eg. "ENSG", "ENTREZ")
#' 
#' @return Results table with gene symbols
#'
#' @examples
#' results = dream_to_dge(model) 
#' results = convert_gene_ids(results, species="human")
#'
#' @importFrom AnnotationDbi mapIds  
#' @importFrom dplyr mutate select rename
#' @export

convert_gene_ids <- function(results, species, id_type = "ENSG") {
  
  # Get appropriate org db for species  
  if(species == "human") {
    org_db <- org.Hs.eg.db
  } else if(species == "mouse") {
    org_db <- org.Mm.eg.db
  }

  if(id_type == "ENSG") {
    # Map Ensembl IDs
    mapped_ids <- mapIds(org_db, keys=results$gene_id,
                         keytype="ENSEMBL", column="SYMBOL") 
  } else if(id_type == "ENTREZ") {
    # Map Entrez IDs
    mapped_ids <- mapIds(org_db, keys=results$gene_id, 
                         keytype="ENTREZID", column="SYMBOL")
  }
  
  # Replace IDs with symbols
  results <- left_join(results, mapped_ids, by="gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::rename(gene_id = SYMBOL)
  
  return(results)
  
}