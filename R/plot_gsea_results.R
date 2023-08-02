#' Visualize GSEA results
#'
#' Generate a forest plot visualization for GSEA results
#'
#' @param gsea.res Data frame of GSEA results 
#' @param padj_cutoff Significance cutoff for shape (default 0.1)
#' @param top_n_pathways Number of up/downregulated pathways to include (default 10)
#' @param remove_first_word_in_pathway_name Remove first word in pathway name (default TRUE)
#'
#' @return ggplot object
#'
#' @examples
#' results <- dge_to_gsea(de_results)
#' plot <- plot_gsea_results(results)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline
#' @importFrom dplyr filter mutate arrange group_by row_number
#' @importFrom stringr str_wrap str_remove word  
#' @importFrom ggforestplot geom_stripes
#' @export

plot_gsea_results <- function(gsea.res, 
                              top_n_pathways = 10, 
                              padj_cutoff = 0.1, 
                              remove_first_word_in_pathway_name = T) {
  
  # Process GSEA results
  processed <- gsea.res %>% 
    clean_msigdbr_pathway_names(remove_first_word = remove_first_word_in_pathway_name) %>% # clean names
    filter(!is.na(NES)) %>%  # remove NA
    group_by(sign(NES)) %>% # group by direction
    arrange(-abs(NES)) %>% # sort by NES
    filter(row_number() <= top_n_pathways) %>% # number of pathways to include
    mutate(
      logP = -log(pval, 10), # p-value
      is_sig = padj < padj_cutoff, # significance
      pathway = str_wrap(pathway, 35) # wrap names
    )
  
  # Generate plot
  plot <- processed %>% 
    ggplot(aes(x = NES, y = reorder(pathway,NES), shape = is_sig)) +
    geom_point(size = 5) + 
    geom_vline(xintercept = 0) +
    ggforestplot::geom_stripes() + # highlight
    scale_shape_manual(values = c("FALSE" = 1,"TRUE" = 16)) +
    theme_classic() +
    theme(aspect.ratio = 3,
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 16)) +
    labs(x = "Net enrichment score", 
         y = "GSEA pathway",
         shape = paste0("padj < ", padj_cutoff)) +
    guides(size = "none")
  
  return(plot)
  
}


# Helper functions

clean_msigdbr_pathway_names <- function(gsea.res, remove_first_word = T){
  
  if(remove_first_word){
    gsea.res <- gsea.res %>%
      mutate(pathway = str_remove(pathway, word(pathway,sep = "_")) %>%
               str_remove("_"))
  }
  
  gsea.res %>%
    mutate(pathway = str_replace_all(pathway, "_", " "))
  
}