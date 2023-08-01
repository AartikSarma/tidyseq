
#' Generate volcano plot for differential expression results
#'
#' @param dge.res Tidy DE results table
#' @param col_up Color for upregulated points (default 'red') 
#' @param col_down Color for downregulated points (default 'blue')
#' @param gene_label_col Column with gene labels to highlight (default NULL) 
#' @param gene_labels Vector of gene labels to highlight (default NULL)  
#' @param n_labels Number of top genes to label (default 10)
#' @param fc_cutoff Fold change cutoff for significance (default 0.5)
#' @param fdr_cutoff FDR cutoff for significance (default 0.1)
#' @param plot_title Plot title text (default '')
#'
#' @return ggplot2 object with volcano plot
#'
#' @examples
#' results <- deseq2_to_dge(dds, contrast=c(1,0))
#' volcano_plot <- plot_volcano(results)
#'
#' @importFrom ggplot2 ggplot aes geom_point 
#' @importFrom dplyr filter arrange case_when mutate
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales squish
#' @export

plot_volcano <- function(dge.res, 
                         col_up = "red", 
                         col_down = "blue",
                         gene_label_col = "hgnc_symbol",
                         gene_labels = NULL,
                         n_labels = 10,
                         fc_cutoff = 0.5,
                         fdr_cutoff = 0.1,
                         plot_title = "") {
  
  dge.res <- 
    dge.res %>%
    mutate(gene_label_col = eval(parse(text = gene_label_col)))
  
  # Label top DE genes if labels not provided
  if(is.null(gene_labels)) {
    gene_labels <- dge.res %>%
      arrange(pvalue) %>%  
      filter(abs(logFC) > fc_cutoff, fdr < fdr_cutoff, !is.na(gene_label_col)) %>%
      group_by(sign(logFC)) %>%
      filter(row_number() <= n_labels) %>%
      pull(gene_label_col)
  }
  
  # Determine default axis limits
  xmin <- min(dge.res$logFC, na.rm = T)
  xmax <- max(dge.res$logFC, na.rm = T)
  ymax <- max(-log(dge.res$pvalue, 10), na.rm = T)
  
  # Color points by significance
  dge.res <- dge.res %>%
    mutate(
      sig.color = case_when(
        fdr < fdr_cutoff & logFC > fc_cutoff ~ col.up,
        fdr < fdr_cutoff & logFC < -fc_cutoff ~ col.down,
        TRUE ~ "gray75"
      ),
      sig.alpha = case_when(
        fdr < fdr_cutoff & logFC > fc_cutoff ~ 1,
        fdr < fdr_cutoff & logFC < -fc_cutoff ~ 1,
        TRUE ~ 0.65
      ),
      labels = case_when(
        (gene_label_col %in% gene_labels) & (fdr < fdr_cutoff) ~ gene_label_col)
    )
  
  # Generate plot
  p <- dge.res %>%
    ggplot(aes(x = logFC, y = -log(pvalue, 10), 
               color = sig.color, alpha = sig.alpha)) +
    geom_point() +
    geom_text_repel(aes(label = labels), color="black") +
    scale_color_identity(guide = "none") + 
    scale_alpha_identity(guide = "none") +
    theme_classic() +
    theme(aspect.ratio = 1, legend.position = "none",
          axis.line = element_blank(),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.border = element_rect(linewidth = 2, fill = NA),
          panel.grid.major = element_line(linetype = "dotted", 
                                          color = "grey85", size = 0.5)) +
    scale_y_continuous(oob = squish, limits = c(0,ymax)) + 
    scale_x_continuous(oob = squish, limits = c(xmin, xmax)) +
    xlab("Log2 fold difference") +
    labs(x=expression(log[2]~fold~difference),  
         y=expression(-log[10]~p~value),
         title = plot.title)
  
  return(p)
  
}