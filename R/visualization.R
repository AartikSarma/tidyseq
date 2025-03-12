#' Generate a PCA plot of RNA-seq data
#'
#' @description Create a PCA plot from normalized expression data
#'
#' @param tidyrna_object A tidyrna object with normalized counts
#' @param color_by Column in metadata to color points by
#' @param shape_by Column in metadata to shape points by (optional)
#' @param label_points Whether to add sample labels to points (default: FALSE)
#' @param components Which principal components to plot (default: c(1, 2))
#' @param interactive Whether to create an interactive plotly plot (default: FALSE)
#'
#' @return A ggplot2 object (or plotly object if interactive = TRUE)
#'
#' @export
plot_pca <- function(tidyrna_object,
                      color_by = NULL,
                      shape_by = NULL,
                      label_points = FALSE,
                      components = c(1, 2),
                      interactive = FALSE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if normalized counts are available
  if (is.null(tidyrna_object$normalized_counts)) {
    stop("No normalized counts found. Run normalize_counts() first.")
  }
  
  # Check if metadata is available
  if (is.null(tidyrna_object$metadata)) {
    stop("No metadata found. Run import_data() first.")
  }
  
  # Log PCA plot creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    "Generating PCA plot",
    "log"
  )
  
  # Get normalized counts and transform if needed
  counts <- tidyrna_object$normalized_counts
  
  # Log transform if not already done
  if (max(counts, na.rm = TRUE) > 30) {
    counts <- log2(counts + 1)
  }
  
  # Transpose to have samples as rows
  counts_t <- t(counts)
  
  # Run PCA
  pca_result <- prcomp(counts_t, scale. = TRUE, center = TRUE)
  
  # Extract PCA data for plotting
  pca_data <- as.data.frame(pca_result$x)
  
  # Add sample information
  pca_data$sample <- rownames(pca_data)
  
  # Add metadata
  metadata <- tidyrna_object$metadata
  pca_data <- merge(pca_data, metadata, by.x = "sample", by.y = "row.names")
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
  
  # Set components to plot
  pc_x <- components[1]
  pc_y <- components[2]
  
  # Create plot labels
  x_lab <- paste0("PC", pc_x, " (", round(var_explained[pc_x], 1), "%)")
  y_lab <- paste0("PC", pc_y, " (", round(var_explained[pc_y], 1), "%)")
  
  # Create base plot
  library(ggplot2)
  
  p <- ggplot(pca_data, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y))) +
    theme_bw() +
    labs(x = x_lab, y = y_lab, title = "PCA Plot") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Add aesthetic mappings
  if (!is.null(color_by)) {
    if (!color_by %in% colnames(pca_data)) {
      stop("Color variable '", color_by, "' not found in metadata")
    }
    p <- p + aes_string(color = color_by)
  }
  
  if (!is.null(shape_by)) {
    if (!shape_by %in% colnames(pca_data)) {
      stop("Shape variable '", shape_by, "' not found in metadata")
    }
    p <- p + aes_string(shape = shape_by)
  }
  
  # Add points
  p <- p + geom_point(size = 3, alpha = 0.8)
  
  # Add sample labels if requested
  if (label_points) {
    p <- p + geom_text(aes(label = sample), hjust = -0.3, vjust = 0.3, size = 3)
  }
  
  # Make interactive if requested
  if (interactive) {
    library(plotly)
    p <- plotly::ggplotly(p)
  }
  
  return(p)
}

#' Generate a heatmap of top differentially expressed genes
#'
#' @description Create a heatmap of top DE genes from a specific contrast
#'
#' @param tidyrna_object A tidyrna object with DE results and normalized counts
#' @param contrast Name of contrast to visualize
#' @param n_genes Number of top genes to include (default: 50)
#' @param cluster_rows Whether to cluster rows (default: TRUE)
#' @param cluster_cols Whether to cluster columns (default: TRUE)
#' @param show_gene_names Whether to show gene names (default: TRUE)
#' @param show_sample_names Whether to show sample names (default: TRUE)
#' @param annotation_cols Columns from metadata to use for annotation (default: NULL)
#' @param color_scale Color scale for heatmap (default: c("blue", "white", "red"))
#'
#' @return A heatmap object
#'
#' @export
plot_de_heatmap <- function(tidyrna_object,
                             contrast,
                             n_genes = 50,
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             show_gene_names = TRUE,
                             show_sample_names = TRUE,
                             annotation_cols = NULL,
                             color_scale = c("blue", "white", "red")) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if DE results are available
  if (is.null(tidyrna_object$de_results) || length(tidyrna_object$de_results) == 0) {
    stop("No DE results found. Run run_de_analysis() first.")
  }
  
  # Check if specified contrast exists
  if (!contrast %in% names(tidyrna_object$de_results)) {
    stop("Specified contrast not found: ", contrast)
  }
  
  # Check if normalized counts are available
  if (is.null(tidyrna_object$normalized_counts)) {
    stop("No normalized counts found. Run normalize_counts() first.")
  }
  
  # Log heatmap creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Generating heatmap for contrast: ", contrast),
    "log"
  )
  
  # Get DE results for the contrast
  de_result <- tidyrna_object$de_results[[contrast]]
  
  # Sort by adjusted p-value
  de_result <- de_result[order(de_result$padj), ]
  
  # Get top genes
  top_genes <- de_result$gene_id[1:min(n_genes, nrow(de_result))]
  
  # Get normalized counts for top genes
  counts <- tidyrna_object$normalized_counts
  top_counts <- counts[rownames(counts) %in% top_genes, ]
  
  # Log transform if not already done
  if (max(top_counts, na.rm = TRUE) > 30) {
    top_counts <- log2(top_counts + 1)
  }
  
  # Scale rows (z-score) for better visualization
  top_counts_scaled <- t(scale(t(top_counts)))
  
  # Create column annotations if requested
  if (!is.null(annotation_cols)) {
    library(pheatmap)
    
    # Check if specified columns exist
    missing_cols <- setdiff(annotation_cols, colnames(tidyrna_object$metadata))
    if (length(missing_cols) > 0) {
      stop("Specified annotation columns not found: ", paste(missing_cols, collapse = ", "))
    }
    
    # Create annotation data frame
    annotation_data <- tidyrna_object$metadata[colnames(top_counts), annotation_cols, drop = FALSE]
    
    # Generate heatmap with annotations
    heatmap <- pheatmap::pheatmap(
      top_counts_scaled,
      annotation_col = annotation_data,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_rownames = show_gene_names,
      show_colnames = show_sample_names,
      color = colorRampPalette(color_scale)(100),
      main = paste0("Top DE Genes: ", contrast),
      fontsize_row = 8,
      fontsize_col = 8
    )
  } else {
    # Generate basic heatmap
    library(pheatmap)
    heatmap <- pheatmap::pheatmap(
      top_counts_scaled,
      cluster_rows = cluster_rows,
      cluster_cols = cluster_cols,
      show_rownames = show_gene_names,
      show_colnames = show_sample_names,
      color = colorRampPalette(color_scale)(100),
      main = paste0("Top DE Genes: ", contrast),
      fontsize_row = 8,
      fontsize_col = 8
    )
  }
  
  return(heatmap)
}

#' Generate a volcano plot for differential expression results
#'
#' @description Create a volcano plot from differential expression analysis
#'
#' @param tidyrna_object A tidyrna object with DE results
#' @param contrast Name of contrast to visualize
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param lfc_cutoff Log fold change cutoff for significance (default: 1)
#' @param label_genes Whether to label top genes (default: TRUE)
#' @param n_labels Number of top genes to label (default: 10)
#' @param interactive Whether to create an interactive plotly plot (default: FALSE)
#'
#' @return A ggplot2 object (or plotly object if interactive = TRUE)
#'
#' @export
plot_volcano <- function(tidyrna_object,
                          contrast,
                          p_cutoff = 0.05,
                          lfc_cutoff = 1,
                          label_genes = TRUE,
                          n_labels = 10,
                          interactive = FALSE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if DE results are available
  if (is.null(tidyrna_object$de_results) || length(tidyrna_object$de_results) == 0) {
    stop("No DE results found. Run run_de_analysis() first.")
  }
  
  # Check if specified contrast exists
  if (!contrast %in% names(tidyrna_object$de_results)) {
    stop("Specified contrast not found: ", contrast)
  }
  
  # Log volcano plot creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Generating volcano plot for contrast: ", contrast),
    "log"
  )
  
  # Get DE results for the contrast
  de_result <- tidyrna_object$de_results[[contrast]]
  
  # Add -log10(pvalue) column
  de_result$neg_log10_pvalue <- -log10(de_result$pvalue)
  
  # Add significance categories
  de_result$category <- "Not significant"
  de_result$category[de_result$padj < p_cutoff & de_result$log2FoldChange > lfc_cutoff] <- "Up-regulated"
  de_result$category[de_result$padj < p_cutoff & de_result$log2FoldChange < -lfc_cutoff] <- "Down-regulated"
  
  # Set factor order for categories
  de_result$category <- factor(
    de_result$category,
    levels = c("Up-regulated", "Down-regulated", "Not significant")
  )
  
  # Create base plot
  library(ggplot2)
  
  p <- ggplot(
    de_result, 
    aes(x = log2FoldChange, y = neg_log10_pvalue, color = category)
  ) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not significant" = "gray")
    ) +
    theme_bw() +
    labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      title = paste0("Volcano Plot: ", contrast),
      color = "Significance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "gray")
  
  # Add gene labels if requested
  if (label_genes) {
    # Find top genes to label (highest combined significance and fold change)
    de_result$score <- de_result$neg_log10_pvalue * abs(de_result$log2FoldChange)
    de_result <- de_result[order(de_result$score, decreasing = TRUE), ]
    top_genes <- de_result[1:min(n_labels, nrow(de_result)), ]
    
    # Add labels to plot
    library(ggrepel)
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene_id),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      min.segment.length = 0.1
    )
  }
  
  # Make interactive if requested
  if (interactive) {
    library(plotly)
    p <- plotly::ggplotly(p)
  }
  
  return(p)
}

#' Generate an enrichment dotplot
#'
#' @description Create a dotplot visualization of enrichment results
#'
#' @param tidyrna_object A tidyrna object with enrichment results
#' @param result_name Name of enrichment result to visualize
#' @param n_terms Number of top terms to include (default: 20)
#' @param color_by Which metric to use for coloring. One of "pvalue", "p.adjust", or "qvalue" (default: "p.adjust")
#' @param size_by Which metric to use for sizing. One of "count" or "GeneRatio" (default: "count")
#' @param interactive Whether to create an interactive plotly plot (default: FALSE)
#'
#' @return A ggplot2 object (or plotly object if interactive = TRUE)
#'
#' @export
plot_enrichment_dotplot <- function(tidyrna_object,
                                     result_name,
                                     n_terms = 20,
                                     color_by = "p.adjust",
                                     size_by = "count",
                                     interactive = FALSE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if enrichment results are available
  if (is.null(tidyrna_object$enrichment_results) || length(tidyrna_object$enrichment_results) == 0) {
    stop("No enrichment results found. Run enrichment functions first.")
  }
  
  # Check if specified result exists
  if (!result_name %in% names(tidyrna_object$enrichment_results)) {
    stop("Specified enrichment result not found: ", result_name)
  }
  
  # Log dotplot creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Generating enrichment dotplot for: ", result_name),
    "log"
  )
  
  # Get enrichment result
  enrichment_result <- tidyrna_object$enrichment_results[[result_name]]
  
  # Handle GO results differently (they're in a list by ontology)
  if (is.list(enrichment_result) && !is.null(enrichment_result$BP)) {
    # This is a GO result with multiple ontologies
    # We'll combine the results and take the top from each ontology
    
    # First, get top terms from each ontology
    top_terms <- list()
    ont_names <- names(enrichment_result)
    
    for (ont in ont_names) {
      if (nrow(enrichment_result[[ont]]@result) > 0) {
        # Sort by adjusted p-value
        result_df <- enrichment_result[[ont]]@result
        result_df <- result_df[order(result_df$p.adjust), ]
        
        # Take top terms
        n_from_ont <- ceiling(n_terms / length(ont_names))
        top_terms[[ont]] <- result_df[1:min(n_from_ont, nrow(result_df)), ]
        top_terms[[ont]]$Ontology <- ont
      }
    }
    
    # Combine results
    if (length(top_terms) > 0) {
      combined_results <- do.call(rbind, top_terms)
      
      # Sort by adjusted p-value and take top terms
      combined_results <- combined_results[order(combined_results$p.adjust), ]
      combined_results <- combined_results[1:min(n_terms, nrow(combined_results)), ]
      
      # Create dotplot
      library(ggplot2)
      
      # Create GeneRatio numeric column
      combined_results$GeneRatio_num <- sapply(
        strsplit(combined_results$GeneRatio, "/"), 
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
      
      # Create base plot
      p <- ggplot(
        combined_results, 
        aes_string(
          x = "Ontology", 
          y = "Description", 
          color = color_by, 
          size = size_by
        )
      ) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        theme_bw() +
        labs(
          x = "GO Ontology",
          y = "Pathway",
          title = paste0("Enrichment Dotplot: ", result_name),
          color = color_by,
          size = size_by
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10)
        )
    } else {
      stop("No enrichment terms found")
    }
  } else {
    # This is a standard enrichment result
    # Make sure we have results
    if (nrow(enrichment_result@result) == 0) {
      stop("No enrichment terms found in result")
    }
    
    # Sort by adjusted p-value
    result_df <- enrichment_result@result
    result_df <- result_df[order(result_df$p.adjust), ]
    
    # Take top terms
    result_df <- result_df[1:min(n_terms, nrow(result_df)), ]
    
    # Create dotplot
    p <- enrichplot::dotplot(enrichment_result, showCategory = n_terms)
  }
  
  # Make interactive if requested
  if (interactive) {
    library(plotly)
    p <- plotly::ggplotly(p)
  }
  
  return(p)
}

#' Generate a MA plot for differential expression results
#'
#' @description Create a MA plot from differential expression analysis
#'
#' @param tidyrna_object A tidyrna object with DE results
#' @param contrast Name of contrast to visualize
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param label_genes Whether to label top genes (default: TRUE)
#' @param n_labels Number of top genes to label (default: 10)
#' @param interactive Whether to create an interactive plotly plot (default: FALSE)
#'
#' @return A ggplot2 object (or plotly object if interactive = TRUE)
#'
#' @export
plot_ma <- function(tidyrna_object,
                     contrast,
                     p_cutoff = 0.05,
                     label_genes = TRUE,
                     n_labels = 10,
                     interactive = FALSE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if DE results are available
  if (is.null(tidyrna_object$de_results) || length(tidyrna_object$de_results) == 0) {
    stop("No DE results found. Run run_de_analysis() first.")
  }
  
  # Check if specified contrast exists
  if (!contrast %in% names(tidyrna_object$de_results)) {
    stop("Specified contrast not found: ", contrast)
  }
  
  # Log MA plot creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Generating MA plot for contrast: ", contrast),
    "log"
  )
  
  # Get DE results for the contrast
  de_result <- tidyrna_object$de_results[[contrast]]
  
  # Calculate the mean expression if not available
  if (!"baseMean" %in% colnames(de_result)) {
    # If we don't have baseMean, calculate it from normalized counts
    if (!is.null(tidyrna_object$normalized_counts)) {
      counts <- tidyrna_object$normalized_counts
      de_result$baseMean <- rowMeans(counts[de_result$gene_id, ])
    } else {
      stop("Cannot create MA plot: baseMean not available in DE results and no normalized counts found")
    }
  }
  
  # Add significance flag
  de_result$significant <- de_result$padj < p_cutoff
  
  # Create base plot
  library(ggplot2)
  
  p <- ggplot(
    de_result, 
    aes(x = log10(baseMean), y = log2FoldChange, color = significant)
  ) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
    ) +
    theme_bw() +
    labs(
      x = "log10(Mean Expression)",
      y = "log2 Fold Change",
      title = paste0("MA Plot: ", contrast),
      color = "Significance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "blue")
  
  # Add gene labels if requested
  if (label_genes) {
    # Find top genes to label (highest combined significance and fold change)
    de_result$score <- -log10(de_result$padj) * abs(de_result$log2FoldChange)
    de_result <- de_result[order(de_result$score, decreasing = TRUE), ]
    de_result <- de_result[!is.na(de_result$score), ]
    top_genes <- de_result[1:min(n_labels, nrow(de_result)), ]
    
    # Add labels to plot
    library(ggrepel)
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene_id),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      min.segment.length = 0.1
    )
  }
  
  # Make interactive if requested
  if (interactive) {
    library(plotly)
    p <- plotly::ggplotly(p)
  }
  
  return(p)
}