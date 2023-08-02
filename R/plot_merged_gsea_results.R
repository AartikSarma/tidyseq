#' Visualize merged GSEA results
#'
#' Create a dot plot to compare GSEA results
#'
#' @param merged.gsea.results Data frame with merged GSEA results 
#' @param padj_cutoff Significance FDR cutoff (default 0.1)
#' @param top_n_pathways Number of top pathways per analysis to include (default 15)
#' @param minimum_significant_groups Number of groups with a significant result for included pathways (default 1 )
#'
#' @return ggplot object
#'
#' @examples
#' merged_results <- merge_gsea(results1, results2)
#' plot <- plot_merged_gsea_results(merged_results)
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradient2 
#' @importFrom dplyr filter mutate arrange group_by ungroup 
#' @importFrom stringr str_wrap  
#' @export

plot_merged_gsea_results <- function(merged_long_gsea_results,
                                     padj_cutoff = 0.1,
                                     top_n_pathways = 15, 
                                     minimum_significant_groups = 1) {
  
  # Process results
  processed <- merged_long_gsea_results %>%
    mutate(is_sig = padj < padj_cutoff) %>% # significance
    filter(!is.na(NES)) %>% # remove NA
    group_by(result_name, sign(NES)) %>% # group
    arrange(-abs(NES)) %>% # sort
    mutate(include_pathway = row_number() <= top_n_pathways) %>% # top n 
    ungroup %>% 
    filter(sum(include_pathway) >= 1, .by = pathway) %>% # keep if top in any
    mutate(NES_prod = prod(NES, na.rm=TRUE), .by="pathway") %>% # NES direction
    clean_msigdbr_pathway_names() %>% # clean names
    mutate(pathway = str_wrap(pathway, 40)) %>%  # wrap names
    filter(sum(is_sig) >= minimum_significant_groups, .by = "pathway") %>%
    mutate(n_sig_groups = sum(is_sig), .by = "pathway")
  # Generate plot
  plot <- processed %>%
    ggplot(aes(x = result_name, y = reorder(pathway, -n_sig_groups), 
               color = NES)) +
    geom_point(aes(shape = is_sig, size = abs(NES))) +
    scale_color_gradient2(high = "red", low = "blue") + 
    scale_size(range = c(2,5)) +
    scale_shape_manual(values = c("TRUE"=16, "FALSE"=21)) +
    theme_classic() +
    theme(aspect.ratio = 4,
          panel.border = element_rect(fill=NA, linewidth=2),
          axis.line = element_blank(),
          axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    labs(x = "", y = "", shape = paste0("padj < ", padj_cutoff)) +
    guides(size = "none")
  
  return(plot)
  
}