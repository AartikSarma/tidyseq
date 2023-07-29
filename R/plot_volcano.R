
#' Generate volcano plot for differential expression results
#'
#' @param deseq.res Tidy DE results table
#' @param col.up Color for upregulated points (default 'red') 
#' @param col.down Color for downregulated points (default 'blue')
#' @param gene.labels Vector of gene labels to highlight (default NULL)  
#' @param n.labels Number of top genes to label (default 10)
#' @param fcCutoff Fold change cutoff for significance (default 0.5)
#' @param fdrCutoff FDR cutoff for significance (default 0.1)
#' @param plot.title Plot title text (default '')
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
#' @export

plot_volcano <- function(deseq.res, 
                         col.up = "red", 
                         col.down = "blue",
                         gene.labels = NULL,
                         n.labels = 10,
                         fcCutoff = 0.5,
                         fdrCutoff = 0.1,
                         plot.title = "") {
  
  # Label top DE genes if labels not provided
  if(is.null(gene.labels)) {
    gene.labels <- deseq.res %>%
      arrange(pvalue) %>%  
      filter(abs(log2FoldChange) > fcCutoff, fdr < fdrCutoff, !is.na(hgnc_symbol)) %>%
      group_by(sign(log2FoldChange)) %>%
      filter(row_number() <= n.labels) %>%
      pull(hgnc_symbol)
  }
  
  # Determine axis limits
  xmin <- min(deseq.res$log2FoldChange, na.rm = T)
  xmax <- max(deseq.res$log2FoldChange, na.rm = T)
  ymax <- max(-log(deseq.res$pvalue, 10), na.rm = T)
  
  # Color points by significance
  deseq.res <- deseq.res %>%
    mutate(
      sig.color = case_when(
        fdr < fdrCutoff & log2FoldChange > fcCutoff ~ col.up,
        fdr < fdrCutoff & log2FoldChange < -fcCutoff ~ col.down,
        TRUE ~ "gray75"
      ),
      sig.alpha = case_when(
        fdr < fdrCutoff & log2FoldChange > fcCutoff ~ 1,
        fdr < fdrCutoff & log2FoldChange < -fcCutoff ~ 1,
        TRUE ~ 0.65
      ),
      labels = case_when(
        (hgnc_symbol %in% gene.labels) & (fdr < fdrCutoff) ~ hgnc_symbol)
    )
  
  # Generate plot
  p <- deseq.res %>%
    ggplot(aes(x = log2FoldChange, y = -log(pvalue, 10), 
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