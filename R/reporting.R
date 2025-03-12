#' Create a quality control report
#'
#' @description Generate a comprehensive quality control report for RNA-Seq data
#'
#' @param tidyrna_object A tidyrna object with raw/filtered/normalized counts
#' @param output_dir Output directory for report
#' @param file_name Report file name (default: "qc_report.html")
#' @param include_sections Sections to include (default: all)
#' @param interactive Whether to include interactive plots (default: TRUE)
#'
#' @return A tidyrna object with QC report path
#'
#' @export
create_qc_report <- function(tidyrna_object,
                              output_dir,
                              file_name = "qc_report.html",
                              include_sections = c("summary", "raw_counts", "filtered_counts", "normalization", "pca", "correlation"),
                              interactive = TRUE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if required data is available
  if (is.null(tidyrna_object$raw_counts)) {
    stop("Raw counts not found. Run import_data() first.")
  }
  
  # Check if output directory exists, create if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create full output path
  output_path <- file.path(output_dir, file_name)
  
  # Log report creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Creating QC report: ", output_path),
    "log"
  )
  
  # Create RMarkdown template
  temp_rmd <- generate_qc_rmd_template(tidyrna_object, include_sections, interactive)
  
  # Define temporary RMarkdown file
  temp_file <- tempfile(fileext = ".Rmd")
  
  # Write template to temporary file
  writeLines(temp_rmd, temp_file)
  
  # Create the report
  rmarkdown::render(
    input = temp_file,
    output_file = file_name,
    output_dir = output_dir,
    quiet = TRUE
  )
  
  # Clean up temporary file
  unlink(temp_file)
  
  # Store report path in tidyrna object
  tidyrna_object$reports$qc_report <- output_path
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("QC report created successfully: ", output_path),
    "log"
  )
  
  return(tidyrna_object)
}

#' Create a differential expression analysis report
#'
#' @description Generate a comprehensive report of differential expression results
#'
#' @param tidyrna_object A tidyrna object with DE results
#' @param contrasts Contrasts to include (default: all)
#' @param output_dir Output directory for report
#' @param file_name Report file name (default: "de_report.html")
#' @param include_sections Sections to include (default: all)
#' @param interactive Whether to include interactive plots (default: TRUE)
#'
#' @return A tidyrna object with DE report path
#'
#' @export
create_de_report <- function(tidyrna_object,
                              contrasts = NULL,
                              output_dir,
                              file_name = "de_report.html",
                              include_sections = c("summary", "volcano", "ma", "heatmap", "top_genes"),
                              interactive = TRUE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if DE results are available
  if (length(tidyrna_object$de_results) == 0) {
    stop("DE results not found. Run run_de_analysis() first.")
  }
  
  # Use all contrasts if not specified
  if (is.null(contrasts)) {
    contrasts <- names(tidyrna_object$de_results)
  } else {
    # Check if specified contrasts exist
    missing_contrasts <- setdiff(contrasts, names(tidyrna_object$de_results))
    if (length(missing_contrasts) > 0) {
      stop("Specified contrasts not found: ", paste(missing_contrasts, collapse = ", "))
    }
  }
  
  # Check if output directory exists, create if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create full output path
  output_path <- file.path(output_dir, file_name)
  
  # Log report creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Creating DE report: ", output_path),
    "log"
  )
  
  # Create RMarkdown template
  temp_rmd <- generate_de_rmd_template(tidyrna_object, contrasts, include_sections, interactive)
  
  # Define temporary RMarkdown file
  temp_file <- tempfile(fileext = ".Rmd")
  
  # Write template to temporary file
  writeLines(temp_rmd, temp_file)
  
  # Create the report
  rmarkdown::render(
    input = temp_file,
    output_file = file_name,
    output_dir = output_dir,
    quiet = TRUE
  )
  
  # Clean up temporary file
  unlink(temp_file)
  
  # Store report path in tidyrna object
  tidyrna_object$reports$de_report <- output_path
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("DE report created successfully: ", output_path),
    "log"
  )
  
  return(tidyrna_object)
}

#' Create an enrichment analysis report
#'
#' @description Generate a comprehensive report of functional enrichment results
#'
#' @param tidyrna_object A tidyrna object with enrichment results
#' @param result_types Enrichment types to include (default: all)
#' @param contrasts Contrasts to include (default: all)
#' @param output_dir Output directory for report
#' @param file_name Report file name (default: "enrichment_report.html")
#' @param interactive Whether to include interactive plots (default: TRUE)
#'
#' @return A tidyrna object with enrichment report path
#'
#' @export
create_enrichment_report <- function(tidyrna_object,
                                      result_types = NULL,
                                      contrasts = NULL,
                                      output_dir,
                                      file_name = "enrichment_report.html",
                                      interactive = TRUE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if enrichment results are available
  if (length(tidyrna_object$enrichment_results) == 0) {
    stop("Enrichment results not found. Run enrichment functions first.")
  }
  
  # Use all result types if not specified
  if (is.null(result_types)) {
    result_types <- sapply(strsplit(names(tidyrna_object$enrichment_results), "_"), function(x) x[2])
    result_types <- unique(result_types)
  }
  
  # Use all contrasts if not specified
  if (is.null(contrasts)) {
    all_names <- names(tidyrna_object$enrichment_results)
    contrast_names <- sapply(strsplit(all_names, "_"), function(x) x[1])
    contrasts <- unique(contrast_names)
  }
  
  # Filter enrichment results to specified types and contrasts
  enrich_results <- list()
  for (name in names(tidyrna_object$enrichment_results)) {
    parts <- strsplit(name, "_")[[1]]
    if (parts[1] %in% contrasts && parts[2] %in% result_types) {
      enrich_results[[name]] <- tidyrna_object$enrichment_results[[name]]
    }
  }
  
  if (length(enrich_results) == 0) {
    stop("No enrichment results found for specified result types and contrasts")
  }
  
  # Check if output directory exists, create if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create full output path
  output_path <- file.path(output_dir, file_name)
  
  # Log report creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Creating enrichment report: ", output_path),
    "log"
  )
  
  # Create RMarkdown template
  temp_rmd <- generate_enrichment_rmd_template(tidyrna_object, enrich_results, interactive)
  
  # Define temporary RMarkdown file
  temp_file <- tempfile(fileext = ".Rmd")
  
  # Write template to temporary file
  writeLines(temp_rmd, temp_file)
  
  # Create the report
  rmarkdown::render(
    input = temp_file,
    output_file = file_name,
    output_dir = output_dir,
    quiet = TRUE
  )
  
  # Clean up temporary file
  unlink(temp_file)
  
  # Store report path in tidyrna object
  tidyrna_object$reports$enrichment_report <- output_path
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Enrichment report created successfully: ", output_path),
    "log"
  )
  
  return(tidyrna_object)
}

#' Create a comprehensive RNA-Seq analysis report
#'
#' @description Generate a comprehensive report including all analysis results
#'
#' @param tidyrna_object A tidyrna object with analysis results
#' @param output_dir Output directory for report
#' @param file_name Report file name (default: "rnaseq_report.html")
#' @param include_sections Sections to include (default: all)
#' @param interactive Whether to include interactive plots (default: TRUE)
#'
#' @return A tidyrna object with full report path
#'
#' @export
create_full_report <- function(tidyrna_object,
                                output_dir,
                                file_name = "rnaseq_report.html",
                                include_sections = c("summary", "qc", "de", "enrichment"),
                                interactive = TRUE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Check if required results exist
  if (is.null(tidyrna_object$raw_counts)) {
    stop("Raw counts not found. Run import_data() first.")
  }
  
  # Check if output directory exists, create if not
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create full output path
  output_path <- file.path(output_dir, file_name)
  
  # Log report creation
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Creating full RNA-Seq report: ", output_path),
    "log"
  )
  
  # Create RMarkdown template
  temp_rmd <- generate_full_rmd_template(tidyrna_object, include_sections, interactive)
  
  # Define temporary RMarkdown file
  temp_file <- tempfile(fileext = ".Rmd")
  
  # Write template to temporary file
  writeLines(temp_rmd, temp_file)
  
  # Create the report
  rmarkdown::render(
    input = temp_file,
    output_file = file_name,
    output_dir = output_dir,
    quiet = TRUE
  )
  
  # Clean up temporary file
  unlink(temp_file)
  
  # Store report path in tidyrna object
  tidyrna_object$reports$full_report <- output_path
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Full RNA-Seq report created successfully: ", output_path),
    "log"
  )
  
  return(tidyrna_object)
}

#' Generate QC report RMarkdown template
#'
#' @param tidyrna_object A tidyrna object
#' @param include_sections Sections to include
#' @param interactive Whether to include interactive plots
#'
#' @return RMarkdown template as character vector
#'
#' @keywords internal
generate_qc_rmd_template <- function(tidyrna_object, include_sections, interactive) {
  
  # Create RMarkdown header
  rmd <- c(
    "---",
    "title: \"RNA-Seq Quality Control Report\"",
    "date: \"`r format(Sys.time(), '%B %d, %Y')`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: cosmo",
    "    highlight: tango",
    "    fig_width: 10",
    "    fig_height: 7",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(dplyr)",
    "library(tidyr)",
    "library(DT)",
    "library(pheatmap)",
    if (interactive) "library(plotly)" else NULL,
    "options(DT.options = list(pageLength = 10, scrollX = TRUE))",
    "```",
    "",
    "# RNA-Seq Quality Control Report",
    "",
    "This report provides a comprehensive quality control analysis of the RNA-Seq data.",
    ""
  )
  
  # Add summary section if included
  if ("summary" %in% include_sections) {
    rmd <- c(rmd,
      "## Data Summary",
      "",
      "```{r summary}",
      "# Get data dimensions",
      "raw_counts <- tidyrna_object$raw_counts",
      "metadata <- tidyrna_object$metadata",
      "has_filtered <- !is.null(tidyrna_object$filtered_counts)",
      "has_normalized <- !is.null(tidyrna_object$normalized_counts)",
      "",
      "# Create summary tables",
      "data_summary <- data.frame(",
      "  Category = c(\"Samples\", \"Genes\", \"Filtered Genes\", \"Normalization Method\"),",
      "  Value = c(",
      "    ncol(raw_counts),",
      "    nrow(raw_counts),",
      "    ifelse(has_filtered, nrow(tidyrna_object$filtered_counts), \"Not performed\"),",
      "    ifelse(has_normalized, tidyrna_object$normalization_method, \"Not performed\")",
      "  )",
      ")",
      "",
      "# Display summary table",
      "DT::datatable(data_summary, options = list(dom = 't'))",
      "```",
      "",
      "### Sample Metadata",
      "",
      "```{r metadata}",
      "DT::datatable(metadata)",
      "```",
      ""
    )
  }
  
  # Add raw counts section if included
  if ("raw_counts" %in% include_sections) {
    rmd <- c(rmd,
      "## Raw Count Distribution",
      "",
      "### Library Sizes",
      "",
      "```{r library_sizes}",
      "# Calculate library sizes",
      "lib_sizes <- colSums(raw_counts)",
      "lib_size_df <- data.frame(Sample = names(lib_sizes), Size = lib_sizes)",
      "",
      "# Create library size plot",
      "p <- ggplot(lib_size_df, aes(x = reorder(Sample, -Size), y = Size / 1e6)) +",
      "  geom_bar(stat = \"identity\", fill = \"steelblue\") +",
      "  theme_bw() +",
      "  labs(x = \"Sample\", y = \"Library Size (millions of reads)\", title = \"Library Sizes\") +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1))",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      "",
      "### Count Distribution",
      "",
      "```{r count_distribution}",
      "# Prepare data for plotting",
      "counts_long <- tidyr::pivot_longer(as.data.frame(raw_counts), ",
      "                                 cols = everything(), ",
      "                                 names_to = \"Sample\", ",
      "                                 values_to = \"Count\")",
      "",
      "# Create boxplot of counts",
      "p <- ggplot(counts_long, aes(x = Sample, y = log10(Count + 1))) +",
      "  geom_boxplot(fill = \"steelblue\") +",
      "  theme_bw() +",
      "  labs(x = \"Sample\", y = \"log10(Count + 1)\", title = \"Raw Count Distribution\") +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1))",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      "",
      "### Density Plot",
      "",
      "```{r density_plot}",
      "# Create density plot of log counts",
      "p <- ggplot(counts_long, aes(x = log10(Count + 1), color = Sample)) +",
      "  geom_density() +",
      "  theme_bw() +",
      "  labs(x = \"log10(Count + 1)\", y = \"Density\", title = \"Raw Count Density\")",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      ""
    )
  }
  
  # Add filtered counts section if included
  if ("filtered_counts" %in% include_sections && !is.null(tidyrna_object$filtered_counts)) {
    rmd <- c(rmd,
      "## Filtered Count Analysis",
      "",
      "```{r filtered_counts}",
      "# Get filtering parameters",
      "filtering_params <- tidyrna_object$parameters$filtering",
      "",
      "# Create parameter info",
      "param_info <- data.frame(",
      "  Parameter = c(\"Minimum Count\", \"Minimum Samples\", \"Count Type\"),",
      "  Value = c(",
      "    filtering_params$min_count,",
      "    filtering_params$min_samples,",
      "    filtering_params$count_type",
      "  )",
      ")",
      "",
      "# Display parameter info",
      "DT::datatable(param_info, options = list(dom = 't'))",
      "```",
      "",
      "### Filtering Results",
      "",
      "```{r filtering_results}",
      "# Calculate filtering statistics",
      "raw_genes <- nrow(raw_counts)",
      "filtered_genes <- nrow(tidyrna_object$filtered_counts)",
      "removed_genes <- raw_genes - filtered_genes",
      "percent_removed <- removed_genes / raw_genes * 100",
      "",
      "# Create results table",
      "filter_results <- data.frame(",
      "  Metric = c(\"Raw Genes\", \"Filtered Genes\", \"Removed Genes\", \"Percent Removed\"),",
      "  Value = c(",
      "    raw_genes,",
      "    filtered_genes,",
      "    removed_genes,",
      "    paste0(round(percent_removed, 2), \"%\")",
      "  )",
      ")",
      "",
      "# Display results table",
      "DT::datatable(filter_results, options = list(dom = 't'))",
      "```",
      "",
      "### Filtered Count Distribution",
      "",
      "```{r filtered_distribution}",
      "# Prepare filtered counts for plotting",
      "filtered_counts <- tidyrna_object$filtered_counts",
      "filtered_long <- tidyr::pivot_longer(as.data.frame(filtered_counts), ",
      "                                   cols = everything(), ",
      "                                   names_to = \"Sample\", ",
      "                                   values_to = \"Count\")",
      "",
      "# Create boxplot of filtered counts",
      "p <- ggplot(filtered_long, aes(x = Sample, y = log10(Count + 1))) +",
      "  geom_boxplot(fill = \"steelblue\") +",
      "  theme_bw() +",
      "  labs(x = \"Sample\", y = \"log10(Count + 1)\", title = \"Filtered Count Distribution\") +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1))",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      ""
    )
  }
  
  # Add normalization section if included
  if ("normalization" %in% include_sections && !is.null(tidyrna_object$normalized_counts)) {
    rmd <- c(rmd,
      "## Normalized Count Analysis",
      "",
      "```{r normalization_info}",
      "# Get normalization parameters",
      "norm_params <- tidyrna_object$parameters$normalization",
      "",
      "# Create parameter info",
      "norm_info <- data.frame(",
      "  Parameter = c(\"Normalization Method\", \"Used Filtered Counts\"),",
      "  Value = c(",
      "    norm_params$method,",
      "    ifelse(norm_params$use_filtered, \"Yes\", \"No\")",
      "  )",
      ")",
      "",
      "# Display parameter info",
      "DT::datatable(norm_info, options = list(dom = 't'))",
      "```",
      "",
      "### Normalized Count Distribution",
      "",
      "```{r normalized_distribution}",
      "# Prepare normalized counts for plotting",
      "norm_counts <- tidyrna_object$normalized_counts",
      "norm_long <- tidyr::pivot_longer(as.data.frame(norm_counts), ",
      "                              cols = everything(), ",
      "                              names_to = \"Sample\", ",
      "                              values_to = \"Count\")",
      "",
      "# Create boxplot of normalized counts",
      "p <- ggplot(norm_long, aes(x = Sample, y = log10(Count + 1))) +",
      "  geom_boxplot(fill = \"steelblue\") +",
      "  theme_bw() +",
      "  labs(x = \"Sample\", y = \"log10(Count + 1)\", title = \"Normalized Count Distribution\") +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1))",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      "",
      "### Normalized Count Density",
      "",
      "```{r normalized_density}",
      "# Create density plot of normalized counts",
      "p <- ggplot(norm_long, aes(x = log10(Count + 1), color = Sample)) +",
      "  geom_density() +",
      "  theme_bw() +",
      "  labs(x = \"log10(Count + 1)\", y = \"Density\", title = \"Normalized Count Density\")",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      ""
    )
  }
  
  # Add PCA section if included
  if ("pca" %in% include_sections && !is.null(tidyrna_object$normalized_counts)) {
    rmd <- c(rmd,
      "## Principal Component Analysis",
      "",
      "```{r pca_setup, echo=FALSE}",
      "# Prepare metadata with colors",
      "metadata <- tidyrna_object$metadata",
      "color_by <- colnames(metadata)[2]  # Use the second column as default color",
      "```",
      "",
      "### PCA Plot (PC1 vs PC2)",
      "",
      "```{r pca_plot}",
      "# Generate PCA plot",
      "pca_plot <- plot_pca(tidyrna_object, color_by = color_by, interactive = ",
      if (interactive) "TRUE" else "FALSE",
      ")",
      "pca_plot",
      "```",
      "",
      "### PCA Plot (PC2 vs PC3)",
      "",
      "```{r pca_plot_23}",
      "# Generate PCA plot for PC2 vs PC3",
      "pca_plot_23 <- plot_pca(tidyrna_object, color_by = color_by, components = c(2, 3), interactive = ",
      if (interactive) "TRUE" else "FALSE",
      ")",
      "pca_plot_23",
      "```",
      ""
    )
  }
  
  # Add correlation section if included
  if ("correlation" %in% include_sections && !is.null(tidyrna_object$normalized_counts)) {
    rmd <- c(rmd,
      "## Sample Correlation Analysis",
      "",
      "```{r correlation}",
      "# Calculate sample correlation matrix",
      "norm_counts <- tidyrna_object$normalized_counts",
      "cor_matrix <- cor(norm_counts, method = \"pearson\")",
      "",
      "# Create heatmap of correlation matrix",
      "# Get metadata for annotations if available",
      "if (!is.null(tidyrna_object$metadata) && ncol(tidyrna_object$metadata) > 1) {",
      "  # Use first two columns of metadata for annotations",
      "  annotation_data <- tidyrna_object$metadata[colnames(cor_matrix), 1:min(2, ncol(tidyrna_object$metadata)), drop = FALSE]",
      "  pheatmap(cor_matrix, main = \"Sample Correlation Heatmap\", annotation_col = annotation_data)",
      "} else {",
      "  # No annotations",
      "  pheatmap(cor_matrix, main = \"Sample Correlation Heatmap\")",
      "}",
      "```",
      ""
    )
  }
  
  # Add footer
  rmd <- c(rmd,
    "---",
    "",
    "## Session Information",
    "",
    "```{r session_info}",
    "sessionInfo()",
    "```"
  )
  
  # Combine into a single string
  rmd_text <- paste(rmd, collapse = "\n")
  
  return(rmd_text)
}

#' Generate DE report RMarkdown template
#'
#' @param tidyrna_object A tidyrna object
#' @param contrasts Contrasts to include
#' @param include_sections Sections to include
#' @param interactive Whether to include interactive plots
#'
#' @return RMarkdown template as character vector
#'
#' @keywords internal
generate_de_rmd_template <- function(tidyrna_object, contrasts, include_sections, interactive) {
  
  # Create RMarkdown header
  rmd <- c(
    "---",
    "title: \"Differential Expression Analysis Report\"",
    "date: \"`r format(Sys.time(), '%B %d, %Y')`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: cosmo",
    "    highlight: tango",
    "    fig_width: 10",
    "    fig_height: 7",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(dplyr)",
    "library(tidyr)",
    "library(DT)",
    "library(pheatmap)",
    if (interactive) "library(plotly)" else NULL,
    "options(DT.options = list(pageLength = 10, scrollX = TRUE))",
    "```",
    "",
    "# Differential Expression Analysis Report",
    "",
    "This report provides a comprehensive analysis of differential expression results.",
    ""
  )
  
  # Add summary section if included
  if ("summary" %in% include_sections) {
    rmd <- c(rmd,
      "## Analysis Summary",
      "",
      "```{r summary}",
      "# Get DE parameters",
      "de_params <- tidyrna_object$parameters$de_analysis",
      "",
      "# Create parameter info",
      "param_info <- data.frame(",
      "  Parameter = c(\"DE Method\", \"Alpha (FDR)\", \"LFC Threshold\"),",
      "  Value = c(",
      "    de_params$method,",
      "    de_params$alpha,",
      "    de_params$lfc_threshold",
      "  )",
      ")",
      "",
      "# Display parameter info",
      "DT::datatable(param_info, options = list(dom = 't'))",
      "```",
      "",
      "### DE Results Summary",
      "",
      "```{r de_summary}",
      "# Create summary table of DE results",
      "contrast_summaries <- list()",
      "",
      "for (contrast in contrasts) {",
      "  de_result <- tidyrna_object$de_results[[contrast]]",
      "  total_genes <- nrow(de_result)",
      "  up_genes <- sum(de_result$log2FoldChange > 0 & de_result$padj < de_params$alpha, na.rm = TRUE)",
      "  down_genes <- sum(de_result$log2FoldChange < 0 & de_result$padj < de_params$alpha, na.rm = TRUE)",
      "  sig_genes <- up_genes + down_genes",
      "  ",
      "  contrast_summaries[[contrast]] <- data.frame(",
      "    Contrast = contrast,",
      "    \"Total Genes\" = total_genes,",
      "    \"Significant Genes\" = sig_genes,",
      "    \"Up-regulated\" = up_genes,",
      "    \"Down-regulated\" = down_genes,",
      "    \"% Significant\" = round(sig_genes / total_genes * 100, 2)",
      "  )",
      "}",
      "",
      "# Combine all summaries",
      "all_summaries <- do.call(rbind, contrast_summaries)",
      "",
      "# Display summary table",
      "DT::datatable(all_summaries)",
      "```",
      ""
    )
  }
  
  # Process each contrast
  for (contrast in contrasts) {
    # Create contrast section
    rmd <- c(rmd,
      paste0("## Contrast: ", contrast),
      ""
    )
    
    # Add volcano plot if included
    if ("volcano" %in% include_sections) {
      rmd <- c(rmd,
        "### Volcano Plot",
        "",
        paste0("```{r volcano_", gsub("[^a-zA-Z0-9]", "_", contrast), "}"),
        "# Generate volcano plot",
        paste0("volcano_plot <- plot_volcano(tidyrna_object, contrast = \"", contrast, "\", interactive = ", 
               if (interactive) "TRUE" else "FALSE", ")"),
        "volcano_plot",
        "```",
        ""
      )
    }
    
    # Add MA plot if included
    if ("ma" %in% include_sections) {
      rmd <- c(rmd,
        "### MA Plot",
        "",
        paste0("```{r ma_", gsub("[^a-zA-Z0-9]", "_", contrast), "}"),
        "# Generate MA plot",
        paste0("ma_plot <- plot_ma(tidyrna_object, contrast = \"", contrast, "\", interactive = ", 
               if (interactive) "TRUE" else "FALSE", ")"),
        "ma_plot",
        "```",
        ""
      )
    }
    
    # Add heatmap if included
    if ("heatmap" %in% include_sections) {
      rmd <- c(rmd,
        "### Heatmap of Top DE Genes",
        "",
        paste0("```{r heatmap_", gsub("[^a-zA-Z0-9]", "_", contrast), "}"),
        "# Generate heatmap",
        paste0("heatmap_plot <- plot_de_heatmap(tidyrna_object, contrast = \"", contrast, "\", n_genes = 50)"),
        "```",
        ""
      )
    }
    
    # Add top genes table if included
    if ("top_genes" %in% include_sections) {
      rmd <- c(rmd,
        "### Top Differentially Expressed Genes",
        "",
        paste0("```{r top_genes_", gsub("[^a-zA-Z0-9]", "_", contrast), "}"),
        "# Get DE results",
        paste0("de_result <- tidyrna_object$de_results[[\\"", contrast, "\\"]]"),
        "",
        "# Sort by adjusted p-value",
        "de_result <- de_result[order(de_result$padj), ]",
        "",
        "# Select columns to display",
        "display_cols <- c(\"gene_id\", \"log2FoldChange\", \"pvalue\", \"padj\", \"baseMean\")",
        "display_cols <- intersect(display_cols, colnames(de_result))",
        "",
        "# Display top genes",
        "DT::datatable(",
        "  de_result[1:min(50, nrow(de_result)), display_cols],",
        "  options = list(scrollX = TRUE)",
        ")",
        "```",
        ""
      )
    }
  }
  
  # Add footer
  rmd <- c(rmd,
    "---",
    "",
    "## Session Information",
    "",
    "```{r session_info}",
    "sessionInfo()",
    "```"
  )
  
  # Combine into a single string
  rmd_text <- paste(rmd, collapse = "\n")
  
  return(rmd_text)
}

#' Generate enrichment report RMarkdown template
#'
#' @param tidyrna_object A tidyrna object
#' @param enrich_results Enrichment results to include
#' @param interactive Whether to include interactive plots
#'
#' @return RMarkdown template as character vector
#'
#' @keywords internal
generate_enrichment_rmd_template <- function(tidyrna_object, enrich_results, interactive) {
  
  # Create RMarkdown header
  rmd <- c(
    "---",
    "title: \"Functional Enrichment Analysis Report\"",
    "date: \"`r format(Sys.time(), '%B %d, %Y')`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: cosmo",
    "    highlight: tango",
    "    fig_width: 10",
    "    fig_height: 7",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(dplyr)",
    "library(tidyr)",
    "library(DT)",
    "library(clusterProfiler)",
    "library(enrichplot)",
    if (interactive) "library(plotly)" else NULL,
    "options(DT.options = list(pageLength = 10, scrollX = TRUE))",
    "```",
    "",
    "# Functional Enrichment Analysis Report",
    "",
    "This report provides a comprehensive analysis of functional enrichment results.",
    ""
  )
  
  # Add summary section
  rmd <- c(rmd,
    "## Analysis Summary",
    "",
    "```{r summary}",
    "# Create summary table of enrichment results",
    "result_summaries <- list()",
    "",
    "for (result_name in names(enrich_results)) {",
    "  result <- enrich_results[[result_name]]",
    "  ",
    "  # Handle GO results differently (they're in a list by ontology)",
    "  if (is.list(result) && !is.null(result$BP)) {",
    "    # This is a GO result with multiple ontologies",
    "    for (ont in names(result)) {",
    "      if (nrow(result[[ont]]@result) > 0) {",
    "        ont_result <- result[[ont]]",
    "        result_summaries[[paste0(result_name, \"_\", ont)]] <- data.frame(",
    "          \"Result Type\" = paste0(result_name, \" (\", ont, \")\"),",
    "          \"Total Terms\" = nrow(ont_result@result),",
    "          \"Significant Terms\" = sum(ont_result@result$p.adjust < 0.05)",
    "        )",
    "      }",
    "    }",
    "  } else {",
    "    # Standard enrichment result",
    "    if (nrow(result@result) > 0) {",
    "      result_summaries[[result_name]] <- data.frame(",
    "        \"Result Type\" = result_name,",
    "        \"Total Terms\" = nrow(result@result),",
    "        \"Significant Terms\" = sum(result@result$p.adjust < 0.05)",
    "      )",
    "    }",
    "  }",
    "}",
    "",
    "# Combine all summaries",
    "if (length(result_summaries) > 0) {",
    "  all_summaries <- do.call(rbind, result_summaries)",
    "  ",
    "  # Display summary table",
    "  DT::datatable(all_summaries)",
    "} else {",
    "  cat(\"No enrichment results available.\")",
    "}",
    "```",
    ""
  )
  
  # Process each enrichment result
  for (result_name in names(enrich_results)) {
    result <- enrich_results[[result_name]]
    
    # Handle GO results differently (they're in a list by ontology)
    if (is.list(result) && !is.null(result$BP)) {
      # This is a GO result with multiple ontologies
      
      # Create result section
      rmd <- c(rmd,
        paste0("## ", result_name),
        ""
      )
      
      # Process each ontology
      for (ont in names(result)) {
        if (nrow(result[[ont]]@result) > 0) {
          # Create ontology section
          rmd <- c(rmd,
            paste0("### ", ont, " Ontology"),
            "",
            "#### Dotplot",
            "",
            paste0("```{r dotplot_", gsub("[^a-zA-Z0-9]", "_", paste0(result_name, "_", ont)), "}"),
            "# Generate dotplot",
            paste0("dotplot(enrich_results[[\\"", result_name, "\\"]][[\\"", ont, "\\"]], showCategory = 20, title = \"", 
                   ont, " Enrichment\")"),
            "```",
            "",
            "#### Enriched Terms",
            "",
            paste0("```{r terms_", gsub("[^a-zA-Z0-9]", "_", paste0(result_name, "_", ont)), "}"),
            "# Get enrichment results",
            paste0("enrichment_result <- enrich_results[[\\"", result_name, "\\"]][[\\"", ont, "\\"]]@result"),
            "",
            "# Display enrichment table",
            "DT::datatable(",
            "  enrichment_result,",
            "  options = list(scrollX = TRUE)",
            ")",
            "```",
            ""
          )
        }
      }
    } else {
      # Standard enrichment result
      if (nrow(result@result) > 0) {
        # Create result section
        rmd <- c(rmd,
          paste0("## ", result_name),
          "",
          "### Dotplot",
          "",
          paste0("```{r dotplot_", gsub("[^a-zA-Z0-9]", "_", result_name), "}"),
          "# Generate dotplot",
          paste0("dotplot(enrich_results[[\\"", result_name, "\\"]], showCategory = 20, title = \"Enrichment\")"),
          "```",
          "",
          "### Enriched Terms",
          "",
          paste0("```{r terms_", gsub("[^a-zA-Z0-9]", "_", result_name), "}"),
          "# Get enrichment results",
          paste0("enrichment_result <- enrich_results[[\\"", result_name, "\\"]]@result"),
          "",
          "# Display enrichment table",
          "DT::datatable(",
          "  enrichment_result,",
          "  options = list(scrollX = TRUE)",
          ")",
          "```",
          ""
        )
      }
    }
  }
  
  # Add footer
  rmd <- c(rmd,
    "---",
    "",
    "## Session Information",
    "",
    "```{r session_info}",
    "sessionInfo()",
    "```"
  )
  
  # Combine into a single string
  rmd_text <- paste(rmd, collapse = "\n")
  
  return(rmd_text)
}

#' Generate comprehensive RNA-Seq report RMarkdown template
#'
#' @param tidyrna_object A tidyrna object
#' @param include_sections Sections to include
#' @param interactive Whether to include interactive plots
#'
#' @return RMarkdown template as character vector
#'
#' @keywords internal
generate_full_rmd_template <- function(tidyrna_object, include_sections, interactive) {
  
  # Create RMarkdown header
  rmd <- c(
    "---",
    "title: \"Comprehensive RNA-Seq Analysis Report\"",
    "date: \"`r format(Sys.time(), '%B %d, %Y')`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: cosmo",
    "    highlight: tango",
    "    fig_width: 10",
    "    fig_height: 7",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(dplyr)",
    "library(tidyr)",
    "library(DT)",
    "library(pheatmap)",
    "library(clusterProfiler)",
    "library(enrichplot)",
    if (interactive) "library(plotly)" else NULL,
    "options(DT.options = list(pageLength = 10, scrollX = TRUE))",
    "```",
    "",
    "# Comprehensive RNA-Seq Analysis Report",
    "",
    "This report provides a comprehensive analysis of RNA-Seq data, including quality control, differential expression, and functional enrichment results.",
    ""
  )
  
  # Add summary section
  rmd <- c(rmd,
    "## Analysis Summary",
    "",
    "```{r summary}",
    "# Get data dimensions",
    "raw_counts <- tidyrna_object$raw_counts",
    "metadata <- tidyrna_object$metadata",
    "has_filtered <- !is.null(tidyrna_object$filtered_counts)",
    "has_normalized <- !is.null(tidyrna_object$normalized_counts)",
    "has_de <- length(tidyrna_object$de_results) > 0",
    "has_enrichment <- length(tidyrna_object$enrichment_results) > 0",
    "",
    "# Create summary tables",
    "data_summary <- data.frame(",
    "  Category = c(\"Samples\", \"Genes\", \"Filtered Genes\", \"Normalization Method\",",
    "               \"DE Contrasts\", \"Enrichment Results\"),",
    "  Value = c(",
    "    ncol(raw_counts),",
    "    nrow(raw_counts),",
    "    ifelse(has_filtered, nrow(tidyrna_object$filtered_counts), \"Not performed\"),",
    "    ifelse(has_normalized, tidyrna_object$normalization_method, \"Not performed\"),",
    "    ifelse(has_de, length(tidyrna_object$de_results), \"Not performed\"),",
    "    ifelse(has_enrichment, length(tidyrna_object$enrichment_results), \"Not performed\")",
    "  )",
    ")",
    "",
    "# Display summary table",
    "DT::datatable(data_summary, options = list(dom = 't'))",
    "```",
    ""
  )
  
  # Add QC section if included
  if ("qc" %in% include_sections) {
    # Include a subset of the QC report
    rmd <- c(rmd,
      "## Quality Control Analysis",
      "",
      "### Sample Overview",
      "",
      "```{r metadata}",
      "# Display metadata",
      "DT::datatable(metadata)",
      "```",
      "",
      "### Library Size Distribution",
      "",
      "```{r library_sizes}",
      "# Calculate library sizes",
      "lib_sizes <- colSums(raw_counts)",
      "lib_size_df <- data.frame(Sample = names(lib_sizes), Size = lib_sizes)",
      "",
      "# Create library size plot",
      "p <- ggplot(lib_size_df, aes(x = reorder(Sample, -Size), y = Size / 1e6)) +",
      "  geom_bar(stat = \"identity\", fill = \"steelblue\") +",
      "  theme_bw() +",
      "  labs(x = \"Sample\", y = \"Library Size (millions of reads)\", title = \"Library Sizes\") +",
      "  theme(axis.text.x = element_text(angle = 90, hjust = 1))",
      "",
      if (interactive) {
        "# Make interactive plot if requested"
        "ggplotly(p)"
      } else {
        "p"
      },
      "```",
      ""
    )
    
    # Add PCA plot if normalized counts are available
    if (!is.null(tidyrna_object$normalized_counts)) {
      rmd <- c(rmd,
        "### Principal Component Analysis",
        "",
        "```{r pca_setup, echo=FALSE}",
        "# Prepare metadata with colors",
        "metadata <- tidyrna_object$metadata",
        "color_by <- colnames(metadata)[2]  # Use the second column as default color",
        "```",
        "",
        "```{r pca_plot}",
        "# Generate PCA plot",
        "pca_plot <- plot_pca(tidyrna_object, color_by = color_by, interactive = ",
        if (interactive) "TRUE" else "FALSE",
        ")",
        "pca_plot",
        "```",
        ""
      )
    }
  }
  
  # Add DE section if included
  if ("de" %in% include_sections && length(tidyrna_object$de_results) > 0) {
    # Get contrasts
    contrasts <- names(tidyrna_object$de_results)
    
    rmd <- c(rmd,
      "## Differential Expression Analysis",
      "",
      "```{r de_summary}",
      "# Get DE parameters",
      "de_params <- tidyrna_object$parameters$de_analysis",
      "",
      "# Create parameter info",
      "param_info <- data.frame(",
      "  Parameter = c(\"DE Method\", \"Alpha (FDR)\", \"LFC Threshold\"),",
      "  Value = c(",
      "    de_params$method,",
      "    de_params$alpha,",
      "    de_params$lfc_threshold",
      "  )",
      ")",
      "",
      "# Display parameter info",
      "DT::datatable(param_info, options = list(dom = 't'))",
      "```",
      "",
      "### DE Results Summary",
      "",
      "```{r de_results_summary}",
      "# Create summary table of DE results",
      "contrast_summaries <- list()",
      "",
      "for (contrast in contrasts) {",
      "  de_result <- tidyrna_object$de_results[[contrast]]",
      "  total_genes <- nrow(de_result)",
      "  up_genes <- sum(de_result$log2FoldChange > 0 & de_result$padj < de_params$alpha, na.rm = TRUE)",
      "  down_genes <- sum(de_result$log2FoldChange < 0 & de_result$padj < de_params$alpha, na.rm = TRUE)",
      "  sig_genes <- up_genes + down_genes",
      "  ",
      "  contrast_summaries[[contrast]] <- data.frame(",
      "    Contrast = contrast,",
      "    \"Total Genes\" = total_genes,",
      "    \"Significant Genes\" = sig_genes,",
      "    \"Up-regulated\" = up_genes,",
      "    \"Down-regulated\" = down_genes,",
      "    \"% Significant\" = round(sig_genes / total_genes * 100, 2)",
      "  )",
      "}",
      "",
      "# Combine all summaries",
      "all_summaries <- do.call(rbind, contrast_summaries)",
      "",
      "# Display summary table",
      "DT::datatable(all_summaries)",
      "```",
      ""
    )
    
    # Add volcano plots for each contrast
    rmd <- c(rmd,
      "### Volcano Plots",
      ""
    )
    
    for (contrast in contrasts) {
      rmd <- c(rmd,
        paste0("#### ", contrast),
        "",
        paste0("```{r volcano_", gsub("[^a-zA-Z0-9]", "_", contrast), "}"),
        "# Generate volcano plot",
        paste0("volcano_plot <- plot_volcano(tidyrna_object, contrast = \"", contrast, "\", interactive = ", 
               if (interactive) "TRUE" else "FALSE", ")"),
        "volcano_plot",
        "```",
        ""
      )
    }
  }
  
  # Add enrichment section if included
  if ("enrichment" %in% include_sections && length(tidyrna_object$enrichment_results) > 0) {
    # Get enrichment results
    enrich_results <- tidyrna_object$enrichment_results
    
    rmd <- c(rmd,
      "## Functional Enrichment Analysis",
      "",
      "```{r enrichment_summary}",
      "# Create summary table of enrichment results",
      "result_summaries <- list()",
      "",
      "for (result_name in names(enrich_results)) {",
      "  result <- enrich_results[[result_name]]",
      "  ",
      "  # Handle GO results differently (they're in a list by ontology)",
      "  if (is.list(result) && !is.null(result$BP)) {",
      "    # This is a GO result with multiple ontologies",
      "    for (ont in names(result)) {",
      "      if (nrow(result[[ont]]@result) > 0) {",
      "        ont_result <- result[[ont]]",
      "        result_summaries[[paste0(result_name, \"_\", ont)]] <- data.frame(",
      "          \"Result Type\" = paste0(result_name, \" (\", ont, \")\"),",
      "          \"Total Terms\" = nrow(ont_result@result),",
      "          \"Significant Terms\" = sum(ont_result@result$p.adjust < 0.05)",
      "        )",
      "      }",
      "    }",
      "  } else {",
      "    # Standard enrichment result",
      "    if (nrow(result@result) > 0) {",
      "      result_summaries[[result_name]] <- data.frame(",
      "        \"Result Type\" = result_name,",
      "        \"Total Terms\" = nrow(result@result),",
      "        \"Significant Terms\" = sum(result@result$p.adjust < 0.05)",
      "      )",
      "    }",
      "  }",
      "}",
      "",
      "# Combine all summaries",
      "if (length(result_summaries) > 0) {",
      "  all_summaries <- do.call(rbind, result_summaries)",
      "  ",
      "  # Display summary table",
      "  DT::datatable(all_summaries)",
      "} else {",
      "  cat(\"No enrichment results available.\")",
      "}",
      "```",
      ""
    )
    
    # Select a subset of enrichment results to display (to avoid making the report too long)
    result_names <- names(enrich_results)
    if (length(result_names) > 4) {
      result_names <- result_names[1:4]
    }
    
    for (result_name in result_names) {
      result <- enrich_results[[result_name]]
      
      # Handle GO results differently (they're in a list by ontology)
      if (is.list(result) && !is.null(result$BP)) {
        # This is a GO result with multiple ontologies
        
        # Create result section
        rmd <- c(rmd,
          paste0("### ", result_name),
          ""
        )
        
        # Process each ontology
        for (ont in names(result)) {
          if (nrow(result[[ont]]@result) > 0) {
            # Create ontology section
            rmd <- c(rmd,
              paste0("#### ", ont, " Ontology"),
              "",
              paste0("```{r dotplot_", gsub("[^a-zA-Z0-9]", "_", paste0(result_name, "_", ont)), "}"),
              "# Generate dotplot",
              paste0("dotplot(enrich_results[[\\"", result_name, "\\"]][[\\"", ont, "\\"]], showCategory = 10, title = \"", 
                     ont, " Enrichment\")"),
              "```",
              ""
            )
          }
        }
      } else {
        # Standard enrichment result
        if (nrow(result@result) > 0) {
          # Create result section
          rmd <- c(rmd,
            paste0("### ", result_name),
            "",
            paste0("```{r dotplot_", gsub("[^a-zA-Z0-9]", "_", result_name), "}"),
            "# Generate dotplot",
            paste0("dotplot(enrich_results[[\\"", result_name, "\\"]], showCategory = 10, title = \"Enrichment\")"),
            "```",
            ""
          )
        }
      }
    }
  }
  
  # Add footer
  rmd <- c(rmd,
    "---",
    "",
    "## Session Information",
    "",
    "```{r session_info}",
    "sessionInfo()",
    "```"
  )
  
  # Combine into a single string
  rmd_text <- paste(rmd, collapse = "\n")
  
  return(rmd_text)
}