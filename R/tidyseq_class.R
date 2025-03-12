#' TidySeq Class
#'
#' @description An S3 class to store RNA-Seq analysis data and results
#'
#' @details The tidyseq object stores all data and results from the RNA-Seq analysis workflow,
#' including count data, metadata, normalization results, differential expression results,
#' and enrichment results.
#'
#' @export
tidyseq <- function() {
  structure(
    list(
      # Data components
      raw_counts = NULL,
      metadata = NULL,
      
      # Preprocessing results
      filtered_counts = NULL,
      normalized_counts = NULL,
      normalization_method = NULL,
      
      # Analysis results
      de_results = list(),
      enrichment_results = list(),
      module_results = list(),
      
      # Parameters and settings
      parameters = list(),
      
      # Reports
      reports = list(),
      
      # Logs and messages
      log = character(),
      messages = character(),
      warnings = character(),
      errors = character()
    ),
    class = "tidyseq"
  )
}

#' Print method for tidyseq objects
#'
#' @param x A tidyseq object
#' @param ... Additional arguments (not used)
#'
#' @export
print.tidyseq <- function(x, ...) {
  cat("TidySeq object\n")
  cat("==============\n\n")
  
  # Data components
  cat("Data:\n")
  if (!is.null(x$raw_counts)) {
    cat(" - Raw counts: ", nrow(x$raw_counts), " genes, ", ncol(x$raw_counts), " samples\n", sep = "")
  } else {
    cat(" - Raw counts: Not loaded\n")
  }
  
  if (!is.null(x$metadata)) {
    cat(" - Metadata: ", nrow(x$metadata), " samples, ", ncol(x$metadata), " variables\n", sep = "")
  } else {
    cat(" - Metadata: Not loaded\n")
  }
  
  # Preprocessing
  cat("\nPreprocessing:\n")
  if (!is.null(x$filtered_counts)) {
    cat(" - Filtered counts: ", nrow(x$filtered_counts), " genes, ", ncol(x$filtered_counts), " samples\n", sep = "")
  } else {
    cat(" - Filtered counts: Not performed\n")
  }
  
  if (!is.null(x$normalized_counts)) {
    cat(" - Normalized counts: ", nrow(x$normalized_counts), " genes, ", 
        ncol(x$normalized_counts), " samples (", x$normalization_method, ")\n", sep = "")
  } else {
    cat(" - Normalized counts: Not performed\n")
  }
  
  # Analysis results
  cat("\nAnalysis Results:\n")
  if (length(x$de_results) > 0) {
    cat(" - Differential expression: ", length(x$de_results), " contrasts\n", sep = "")
  } else {
    cat(" - Differential expression: Not performed\n")
  }
  
  if (length(x$enrichment_results) > 0) {
    cat(" - Enrichment analysis: ", length(x$enrichment_results), " results\n", sep = "")
  } else {
    cat(" - Enrichment analysis: Not performed\n")
  }
  
  if (length(x$module_results) > 0) {
    cat(" - Module analysis: ", length(x$module_results), " results\n", sep = "")
  } else {
    cat(" - Module analysis: Not performed\n")
  }
  
  # Reports
  cat("\nReports:\n")
  if (length(x$reports) > 0) {
    for (report_name in names(x$reports)) {
      cat(" - ", report_name, ": ", x$reports[[report_name]], "\n", sep = "")
    }
  } else {
    cat(" - No reports generated\n")
  }
}

#' Summary method for tidyseq objects
#'
#' @param object A tidyseq object
#' @param ... Additional arguments (not used)
#'
#' @export
summary.tidyseq <- function(object, ...) {
  print(object)
  
  # Add additional summary information here as needed
  if (length(object$de_results) > 0) {
    cat("\nDifferential Expression Summary:\n")
    for (contrast_name in names(object$de_results)) {
      de_result <- object$de_results[[contrast_name]]
      up_genes <- sum(de_result$log2FoldChange > 0 & de_result$padj < 0.05, na.rm = TRUE)
      down_genes <- sum(de_result$log2FoldChange < 0 & de_result$padj < 0.05, na.rm = TRUE)
      cat(" - ", contrast_name, ": ", up_genes, " up-regulated, ", down_genes, 
          " down-regulated (padj < 0.05)\n", sep = "")
    }
  }
  
  if (length(object$warnings) > 0) {
    cat("\nWarnings:\n")
    for (warning in object$warnings) {
      cat(" - ", warning, "\n", sep = "")
    }
  }
}

#' Get method for tidyseq objects
#'
#' @param tidyseq_object A tidyseq object
#' @param component The component to retrieve
#'
#' @return The requested component from the tidyseq object
#'
#' @export
get_component <- function(tidyseq_object, component) {
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  if (!component %in% names(tidyseq_object)) {
    stop("Component '", component, "' not found in tidyseq object")
  }
  
  return(tidyseq_object[[component]])
}

#' Add log message to tidyseq object
#'
#' @param tidyseq_object A tidyseq object
#' @param message The message to add to the log
#' @param type The type of message (log, message, warning, error)
#'
#' @return Updated tidyseq object with added message
#'
#' @export
add_message <- function(tidyseq_object, message, type = "log") {
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  if (!type %in% c("log", "message", "warning", "error")) {
    stop("Message type must be one of: log, message, warning, error")
  }
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- paste0("[", timestamp, "] ", message)
  
  tidyseq_object[[type]] <- c(tidyseq_object[[type]], formatted_message)
  
  if (type %in% c("warning", "error")) {
    if (type == "warning") {
      warning(message)
    } else {
      stop(message)
    }
  }
  
  return(tidyseq_object)
}