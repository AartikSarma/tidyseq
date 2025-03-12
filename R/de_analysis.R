#' Create contrast specifications for differential expression analysis
#'
#' @description Create contrast specifications from conditions or manually specified contrasts
#'
#' @param tidyrna_object A tidyrna object with metadata
#' @param condition_column Column in metadata containing condition information
#' @param reference_level Reference level for condition comparisons
#' @param contrasts List of manual contrasts as character vectors ("Group1-Group2")
#' @param coef Coefficient numbers/names from design matrix to test (instead of contrasts)
#' @param design_formula Design formula for the model (default: ~ condition)
#'
#' @return A tidyrna object with contrast specifications
#'
#' @export
create_de_contrasts <- function(tidyrna_object,
                                 condition_column = NULL,
                                 reference_level = NULL,
                                 contrasts = NULL,
                                 coef = NULL,
                                 design_formula = NULL) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  if (is.null(tidyrna_object$metadata)) {
    stop("No metadata found in tidyrna object. Run import_data() first.")
  }
  
  # Initialize contrasts list
  contrast_list <- list()
  
  # Create design formula if not provided
  if (is.null(design_formula) && !is.null(condition_column)) {
    design_formula <- as.formula(paste0("~ ", condition_column))
  }
  
  # Case 1: Use condition column and reference level
  if (!is.null(condition_column)) {
    tidyrna_object <- add_message(
      tidyrna_object, 
      paste0("Creating contrasts from condition column: ", condition_column),
      "log"
    )
    
    # Check if condition column exists
    if (!condition_column %in% colnames(tidyrna_object$metadata)) {
      stop("Condition column '", condition_column, "' not found in metadata")
    }
    
    # Get condition levels
    conditions <- as.factor(tidyrna_object$metadata[[condition_column]])
    levels <- levels(conditions)
    
    # Set reference level if specified
    if (!is.null(reference_level)) {
      if (!reference_level %in% levels) {
        stop("Reference level '", reference_level, "' not found in condition levels")
      }
      
      # Reorder levels to make reference level first
      levels <- c(reference_level, levels[levels != reference_level])
      conditions <- factor(conditions, levels = levels)
    }
    
    # Create all pairwise contrasts
    if (is.null(contrasts)) {
      contrasts <- list()
      for (i in 2:length(levels)) {
        contrast_name <- paste0(levels[i], "_vs_", levels[1])
        contrasts[[contrast_name]] <- c(levels[i], levels[1])
      }
    }
    
    # Store design formula and condition information
    tidyrna_object$parameters$de_design <- list(
      formula = design_formula,
      condition_column = condition_column,
      reference_level = reference_level,
      levels = levels
    )
    
    contrast_list <- contrasts
  }
  
  # Case 2: Use manually specified contrasts
  else if (!is.null(contrasts)) {
    tidyrna_object <- add_message(
      tidyrna_object, 
      "Using manually specified contrasts",
      "log"
    )
    
    contrast_list <- contrasts
  }
  
  # Case 3: Use coefficients
  else if (!is.null(coef)) {
    tidyrna_object <- add_message(
      tidyrna_object, 
      "Using specified coefficients for DE analysis",
      "log"
    )
    
    # Store coefficient information
    tidyrna_object$parameters$de_design <- list(
      formula = design_formula,
      coef = coef
    )
  }
  
  # Store contrasts in object
  tidyrna_object$parameters$contrasts <- contrast_list
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    paste0("Created ", length(contrast_list), " contrast specifications"),
    "log"
  )
  
  return(tidyrna_object)
}

#' Run differential expression analysis
#'
#' @description Perform differential expression analysis using DESeq2 or limma
#'
#' @param tidyrna_object A tidyrna object with count data and metadata
#' @param method Analysis method. One of "DESeq2", "limma-voom", or "edgeR" (default: "DESeq2")
#' @param design_formula Design formula (default: from create_de_contrasts)
#' @param contrasts List of contrasts to test (default: from create_de_contrasts)
#' @param coef Coefficient numbers/names to test (default: from create_de_contrasts)
#' @param alpha Significance level cutoff (default: 0.05)
#' @param lfc_threshold Log fold change threshold (default: 0)
#' @param use_filtered Whether to use filtered counts (default: TRUE)
#' @param use_normalized Whether to use pre-computed normalized counts (default: TRUE)
#'
#' @return A tidyrna object with DE results
#'
#' @export
run_de_analysis <- function(tidyrna_object,
                             method = "DESeq2",
                             design_formula = NULL,
                             contrasts = NULL,
                             coef = NULL,
                             alpha = 0.05,
                             lfc_threshold = 0,
                             use_filtered = TRUE,
                             use_normalized = TRUE) {
  
  if (!inherits(tidyrna_object, "tidyrna")) {
    stop("Object must be of class 'tidyrna'")
  }
  
  # Get design information
  if (is.null(design_formula) && !is.null(tidyrna_object$parameters$de_design)) {
    design_formula <- tidyrna_object$parameters$de_design$formula
  }
  
  if (is.null(contrasts) && !is.null(tidyrna_object$parameters$contrasts)) {
    contrasts <- tidyrna_object$parameters$contrasts
  }
  
  if (is.null(coef) && !is.null(tidyrna_object$parameters$de_design$coef)) {
    coef <- tidyrna_object$parameters$de_design$coef
  }
  
  # Check if we have design information
  if (is.null(design_formula) && is.null(contrasts) && is.null(coef)) {
    stop("No design formula, contrasts, or coefficients specified. Run create_de_contrasts() first.")
  }
  
  # Check if we have counts
  if (use_filtered && is.null(tidyrna_object$filtered_counts)) {
    stop("No filtered counts found. Run filter_low_counts() first.")
  }
  
  if (!use_filtered && is.null(tidyrna_object$raw_counts)) {
    stop("No raw counts found. Run import_data() first.")
  }
  
  # Get appropriate count data
  counts <- if (use_filtered) tidyrna_object$filtered_counts else tidyrna_object$raw_counts
  
  # Get metadata
  if (is.null(tidyrna_object$metadata)) {
    stop("No metadata found. Run import_data() first.")
  }
  
  metadata <- tidyrna_object$metadata
  
  # Ensure samples match
  common_samples <- intersect(colnames(counts), rownames(metadata))
  if (length(common_samples) == 0) {
    stop("No common samples found between counts and metadata")
  }
  
  counts <- counts[, common_samples]
  metadata <- metadata[common_samples, ]
  
  # Run DE analysis based on method
  if (method == "DESeq2") {
    tidyrna_object <- run_deseq2_analysis(
      tidyrna_object, counts, metadata, design_formula, 
      contrasts, coef, alpha, lfc_threshold, use_normalized
    )
  } else if (method == "limma-voom") {
    tidyrna_object <- run_limma_voom_analysis(
      tidyrna_object, counts, metadata, design_formula, 
      contrasts, coef, alpha, lfc_threshold, use_normalized
    )
  } else if (method == "edgeR") {
    tidyrna_object <- run_edger_analysis(
      tidyrna_object, counts, metadata, design_formula, 
      contrasts, coef, alpha, lfc_threshold, use_normalized
    )
  } else {
    stop("Unsupported DE analysis method: ", method)
  }
  
  # Store DE analysis parameters
  tidyrna_object$parameters$de_analysis <- list(
    method = method,
    alpha = alpha,
    lfc_threshold = lfc_threshold,
    use_filtered = use_filtered,
    use_normalized = use_normalized
  )
  
  return(tidyrna_object)
}

#' Run DESeq2 differential expression analysis
#'
#' @param tidyrna_object A tidyrna object
#' @param counts Count data
#' @param metadata Sample metadata
#' @param design_formula Design formula
#' @param contrasts List of contrasts to test
#' @param coef Coefficient numbers/names to test
#' @param alpha Significance level cutoff
#' @param lfc_threshold Log fold change threshold
#' @param use_normalized Whether to use pre-computed normalized counts
#'
#' @return Updated tidyrna object with DESeq2 results
#'
#' @keywords internal
run_deseq2_analysis <- function(tidyrna_object, 
                                counts, 
                                metadata, 
                                design_formula, 
                                contrasts, 
                                coef, 
                                alpha, 
                                lfc_threshold,
                                use_normalized) {
  
  library(DESeq2)
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    "Running DESeq2 differential expression analysis",
    "log"
  )
  
  # Use existing DESeq2 object if available and appropriate
  if (use_normalized && !is.null(tidyrna_object$deseq2_object)) {
    dds <- tidyrna_object$deseq2_object
    
    # Update design if necessary
    if (!is.null(design_formula) && design_formula != design(dds)) {
      design(dds) <- design_formula
      dds <- DESeq2::estimateDispersions(dds)
    }
    
  } else {
    # Create new DESeq2 object
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = metadata,
      design = design_formula
    )
    
    # Run DESeq2 analysis
    dds <- DESeq2::DESeq(dds)
  }
  
  # Store results for each contrast or coefficient
  results_list <- list()
  
  if (!is.null(contrasts)) {
    # Process each contrast
    for (contrast_name in names(contrasts)) {
      contrast_pair <- contrasts[[contrast_name]]
      
      # Extract results
      res <- DESeq2::results(
        dds, 
        contrast = c(gsub("~", "", deparse(design_formula)[[1]]), contrast_pair[1], contrast_pair[2]),
        alpha = alpha,
        lfcThreshold = lfc_threshold
      )
      
      # Convert to data frame
      res_df <- as.data.frame(res)
      res_df$gene_id <- rownames(res_df)
      
      # Add significance flag
      res_df$significant <- res_df$padj < alpha & abs(res_df$log2FoldChange) > lfc_threshold
      
      # Store results
      results_list[[contrast_name]] <- res_df
      
      # Log results
      sig_count <- sum(res_df$significant, na.rm = TRUE)
      tidyrna_object <- add_message(
        tidyrna_object, 
        paste0("Contrast ", contrast_name, ": ", sig_count, " significant genes (padj < ", 
               alpha, ", |log2FC| > ", lfc_threshold, ")"),
        "log"
      )
    }
  } else if (!is.null(coef)) {
    # Process each coefficient
    for (coef_name in coef) {
      # Extract results
      res <- DESeq2::results(
        dds, 
        name = coef_name,
        alpha = alpha,
        lfcThreshold = lfc_threshold
      )
      
      # Convert to data frame
      res_df <- as.data.frame(res)
      res_df$gene_id <- rownames(res_df)
      
      # Add significance flag
      res_df$significant <- res_df$padj < alpha & abs(res_df$log2FoldChange) > lfc_threshold
      
      # Store results
      results_list[[coef_name]] <- res_df
      
      # Log results
      sig_count <- sum(res_df$significant, na.rm = TRUE)
      tidyrna_object <- add_message(
        tidyrna_object, 
        paste0("Coefficient ", coef_name, ": ", sig_count, " significant genes (padj < ", 
               alpha, ", |log2FC| > ", lfc_threshold, ")"),
        "log"
      )
    }
  } else {
    # No contrasts or coefficients specified
    warning("No contrasts or coefficients specified. No DE results generated.")
  }
  
  # Store results in tidyrna object
  tidyrna_object$de_results <- results_list
  
  # Store DESeq2 object
  tidyrna_object$deseq2_object <- dds
  
  return(tidyrna_object)
}

#' Run limma-voom differential expression analysis
#'
#' @param tidyrna_object A tidyrna object
#' @param counts Count data
#' @param metadata Sample metadata
#' @param design_formula Design formula
#' @param contrasts List of contrasts to test
#' @param coef Coefficient numbers/names to test
#' @param alpha Significance level cutoff
#' @param lfc_threshold Log fold change threshold
#' @param use_normalized Whether to use pre-computed normalized counts
#'
#' @return Updated tidyrna object with limma-voom results
#'
#' @keywords internal
run_limma_voom_analysis <- function(tidyrna_object, 
                                    counts, 
                                    metadata, 
                                    design_formula, 
                                    contrasts, 
                                    coef, 
                                    alpha, 
                                    lfc_threshold,
                                    use_normalized) {
  
  library(limma)
  library(edgeR)
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    "Running limma-voom differential expression analysis",
    "log"
  )
  
  # Create design matrix
  design <- model.matrix(design_formula, data = metadata)
  
  # Create DGEList object
  dge <- edgeR::DGEList(counts = counts)
  
  # Calculate normalization factors
  dge <- edgeR::calcNormFactors(dge)
  
  # Run voom
  v <- limma::voom(dge, design, plot = FALSE)
  
  # Fit linear model
  fit <- limma::lmFit(v, design)
  
  # Store results for each contrast or coefficient
  results_list <- list()
  
  if (!is.null(contrasts)) {
    # Process each contrast
    contrast_matrix <- limma::makeContrasts(contrasts = contrasts, levels = design)
    
    # Fit contrasts
    contrast_fit <- limma::contrasts.fit(fit, contrast_matrix)
    contrast_fit <- limma::eBayes(contrast_fit)
    
    # Extract results for each contrast
    for (contrast_name in colnames(contrast_matrix)) {
      # Extract results
      res <- limma::topTable(
        contrast_fit, 
        coef = contrast_name,
        number = Inf,
        sort.by = "P"
      )
      
      # Add gene ID
      res$gene_id <- rownames(res)
      
      # Rename columns to match DESeq2 output
      colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
      colnames(res)[colnames(res) == "P.Value"] <- "pvalue"
      colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"
      
      # Add significance flag
      res$significant <- res$padj < alpha & abs(res$log2FoldChange) > lfc_threshold
      
      # Store results
      results_list[[contrast_name]] <- res
      
      # Log results
      sig_count <- sum(res$significant, na.rm = TRUE)
      tidyrna_object <- add_message(
        tidyrna_object, 
        paste0("Contrast ", contrast_name, ": ", sig_count, " significant genes (padj < ", 
               alpha, ", |log2FC| > ", lfc_threshold, ")"),
        "log"
      )
    }
  } else if (!is.null(coef)) {
    # Fit eBayes
    fit <- limma::eBayes(fit)
    
    # Process each coefficient
    for (coef_name in coef) {
      # Extract results
      res <- limma::topTable(
        fit, 
        coef = coef_name,
        number = Inf,
        sort.by = "P"
      )
      
      # Add gene ID
      res$gene_id <- rownames(res)
      
      # Rename columns to match DESeq2 output
      colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
      colnames(res)[colnames(res) == "P.Value"] <- "pvalue"
      colnames(res)[colnames(res) == "adj.P.Val"] <- "padj"
      
      # Add significance flag
      res$significant <- res$padj < alpha & abs(res$log2FoldChange) > lfc_threshold
      
      # Store results
      results_list[[coef_name]] <- res
      
      # Log results
      sig_count <- sum(res$significant, na.rm = TRUE)
      tidyrna_object <- add_message(
        tidyrna_object, 
        paste0("Coefficient ", coef_name, ": ", sig_count, " significant genes (padj < ", 
               alpha, ", |log2FC| > ", lfc_threshold, ")"),
        "log"
      )
    }
  } else {
    # No contrasts or coefficients specified
    warning("No contrasts or coefficients specified. No DE results generated.")
  }
  
  # Store results in tidyrna object
  tidyrna_object$de_results <- results_list
  
  # Store limma-voom objects
  tidyrna_object$limma_voom_object <- list(
    dge = dge,
    voom = v,
    fit = fit
  )
  
  return(tidyrna_object)
}

#' Run edgeR differential expression analysis
#'
#' @param tidyrna_object A tidyrna object
#' @param counts Count data
#' @param metadata Sample metadata
#' @param design_formula Design formula
#' @param contrasts List of contrasts to test
#' @param coef Coefficient numbers/names to test
#' @param alpha Significance level cutoff
#' @param lfc_threshold Log fold change threshold
#' @param use_normalized Whether to use pre-computed normalized counts
#'
#' @return Updated tidyrna object with edgeR results
#'
#' @keywords internal
run_edger_analysis <- function(tidyrna_object, 
                               counts, 
                               metadata, 
                               design_formula, 
                               contrasts, 
                               coef, 
                               alpha, 
                               lfc_threshold,
                               use_normalized) {
  
  library(edgeR)
  
  tidyrna_object <- add_message(
    tidyrna_object, 
    "Running edgeR differential expression analysis",
    "log"
  )
  
  # Use existing edgeR object if available and appropriate
  if (use_normalized && !is.null(tidyrna_object$edger_object)) {
    dge <- tidyrna_object$edger_object
  } else {
    # Create DGEList object
    dge <- edgeR::DGEList(counts = counts)
    
    # Calculate normalization factors
    dge <- edgeR::calcNormFactors(dge)
  }
  
  # Create design matrix
  design <- model.matrix(design_formula, data = metadata)
  
  # Estimate dispersion
  dge <- edgeR::estimateDisp(dge, design)
  
  # Fit model
  fit <- edgeR::glmQLFit(dge, design)
  
  # Store results for each contrast or coefficient
  results_list <- list()
  
  if (!is.null(contrasts)) {
    # Process each contrast
    for (contrast_name in names(contrasts)) {
      contrast_pair <- contrasts[[contrast_name]]
      
      # Create contrast
      col_names <- gsub("^.*\\((.*)\\)$", "\\1", colnames(design))
      col_idx1 <- which(col_names == contrast_pair[1])
      col_idx2 <- which(col_names == contrast_pair[2])
      
      if (length(col_idx1) == 0 || length(col_idx2) == 0) {
        warning("Could not find contrast levels in design matrix: ", 
                paste(contrast_pair, collapse = " vs "))
        next
      }
      
      contrast <- rep(0, ncol(design))
      contrast[col_idx1] <- 1
      contrast[col_idx2] <- -1
      
      # Test contrast
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
      
      # Extract results
      res <- edgeR::topTags(qlf, n = Inf)$table
      
      # Add gene ID
      res$gene_id <- rownames(res)
      
      # Rename columns to match DESeq2 output
      colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
      colnames(res)[colnames(res) == "PValue"] <- "pvalue"
      colnames(res)[colnames(res) == "FDR"] <- "padj"
      
      # Add significance flag
      res$significant <- res$padj < alpha & abs(res$log2FoldChange) > lfc_threshold
      
      # Store results
      results_list[[contrast_name]] <- res
      
      # Log results
      sig_count <- sum(res$significant, na.rm = TRUE)
      tidyrna_object <- add_message(
        tidyrna_object, 
        paste0("Contrast ", contrast_name, ": ", sig_count, " significant genes (padj < ", 
               alpha, ", |log2FC| > ", lfc_threshold, ")"),
        "log"
      )
    }
  } else if (!is.null(coef)) {
    # Process each coefficient
    for (coef_name in coef) {
      # Get coefficient index
      coef_idx <- if (is.numeric(coef_name)) coef_name else which(colnames(design) == coef_name)
      
      if (length(coef_idx) == 0) {
        warning("Could not find coefficient in design matrix: ", coef_name)
        next
      }
      
      # Test coefficient
      qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)
      
      # Extract results
      res <- edgeR::topTags(qlf, n = Inf)$table
      
      # Add gene ID
      res$gene_id <- rownames(res)
      
      # Rename columns to match DESeq2 output
      colnames(res)[colnames(res) == "logFC"] <- "log2FoldChange"
      colnames(res)[colnames(res) == "PValue"] <- "pvalue"
      colnames(res)[colnames(res) == "FDR"] <- "padj"
      
      # Add significance flag
      res$significant <- res$padj < alpha & abs(res$log2FoldChange) > lfc_threshold
      
      # Store results
      results_list[[as.character(coef_name)]] <- res
      
      # Log results
      sig_count <- sum(res$significant, na.rm = TRUE)
      tidyrna_object <- add_message(
        tidyrna_object, 
        paste0("Coefficient ", coef_name, ": ", sig_count, " significant genes (padj < ", 
               alpha, ", |log2FC| > ", lfc_threshold, ")"),
        "log"
      )
    }
  } else {
    # No contrasts or coefficients specified
    warning("No contrasts or coefficients specified. No DE results generated.")
  }
  
  # Store results in tidyrna object
  tidyrna_object$de_results <- results_list
  
  # Store edgeR objects
  tidyrna_object$edger_object <- list(
    dge = dge,
    fit = fit
  )
  
  return(tidyrna_object)
}