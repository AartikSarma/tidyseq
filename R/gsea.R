#' Perform gene set enrichment analysis
#'
#' Run GSEA on differential expression results using msigdb gene sets
#'
#' @param dge.res Tidy DE results table with fold changes
#' @param nPermSimple Number of permutations for p-value calc (default 10000)  
#' @param sig.list Optional named list of gene sets to use (default NULL)
#' @param species Species for msigdb sets (default 'Homo sapiens')
#' @param msigdb_category msigdb category to use (default 'H')
#' @param msigdb_subcategory msigdb subcategory (default NULL)
#' @param minSize Minimum gene set size (default 30)
#' @param ... Additional parameters for fgseaMultilevel()
#'
#' @return DataFrame with GSEA results
#'
#' @examples
#' results <- deseq2_to_dge(dds)
#' gsea_res <- dge_to_gsea(results)
#'
#' @importFrom fgsea fgseaMultilevel 
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr arrange filter select summarize group_by
#' @importFrom tibble deframe
#' @export

dge_to_gsea <- function(dge.res,
                        nPermSimple = 10000, 
                        sig.list = NULL,
                        species = "Homo sapiens",
                        msigdb_category = "H",
                        msigdb_subcategory = NULL,
                        minSize = 30,
                        ...) {
  
  # Get msigdb gene sets
  if(is.null(sig.list)) {  
    all_gene_sets <- msigdbr::msigdbr(species = species,
                                      category = msigdb_category,
                                      subcategory = msigdb_subcategory)
    msigdbr_list <- split(x = all_gene_sets$gene_symbol, 
                          f = all_gene_sets$gs_name)
  } else {
    msigdbr_list <- sig.list 
  }
  
  # Prepare DE results
  de_mat <- dge.res %>%
    arrange(-logFC) %>%
    filter(!is.na(logFC), !is.na(hgnc_symbol)) %>%  
    select(hgnc_symbol, logFC) %>%
    group_by(hgnc_symbol) %>%
    summarize(logFC = mean(logFC)) %>%
    tibble::deframe()
  
  # Run GSEA 
  gsea_res <- fgsea::fgseaMultilevel(pathways = msigdbr_list, de_mat,
                                     minSize = minSize, 
                                     nPermSimple = nPermSimple,
                                     BPPARAM = BiocParallel::SnowParam(6, progressbar = T),
                                     ...)
  
  return(gsea_res)
  
}