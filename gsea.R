gsea <- function(res,
                 nPermSimple = 10000,
                 sig.list = NULL,
                 species = "Homo sapiens", 
                 msigdb_category = "H", 
                 msigdb_subcategory = NULL, 
                 minSize = 30,
                 ...){
  
  if(is.null(sig.list)){
    all_gene_sets = msigdbr(
      species = species, 
      category = msigdb_category,
      subcategory =  msigdb_subcategory)
    
    msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)
  }
  
  if(!is.null(sig.list)){
    msigdbr_list = sig.list
  }
  
  res %>% 
    arrange(-log2FoldChange) %>%
    filter(!is.na(log2FoldChange), !is.na(hgnc_symbol)) %>%
    dplyr::select(hgnc_symbol, log2FoldChange) %>%
    group_by(hgnc_symbol) %>%
    summarize(log2FoldChange = mean(log2FoldChange) ) %>%
    deframe %>%
    fgseaMultilevel(pathways = msigdbr_list, ., 
                    minSize = minSize,nPermSimple = nPermSimple,
                    BPPARAM =  SnowParam(6, progressbar = T), ...) 
}