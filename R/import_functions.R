#' Import count data and metadata for RNA-Seq analysis
#'
#' @description Import RNA-Seq count data and metadata for analysis with TidyRNA
#'
#' @param counts_file Path to count data file or a data frame/matrix of counts
#' @param metadata_file Path to metadata file or a data frame of metadata
#' @param counts_format Format of the counts file. One of "csv", "tsv", "txt", or "auto"
#' @param metadata_format Format of the metadata file. One of "csv", "tsv", "txt", or "auto"
#' @param sample_column Column in metadata that identifies sample IDs (default: first column)
#' @param gene_column Column in counts data that identifies gene IDs (default: first column)
#' @param validate_ids Whether to validate matching of sample IDs (default: TRUE)
#'
#' @return A tidyseq object with imported data
#'
#' @export
import_data <- function(counts_file, 
                         metadata_file,
                         counts_format = "auto",
                         metadata_format = "auto",
                         sample_column = NULL,
                         gene_column = NULL,
                         validate_ids = TRUE) {
  
  # Initialize tidyseq object
  tidyseq_obj <- tidyseq()
  
  # Import counts
  tidyseq_obj <- add_message(tidyseq_obj, "Importing count data", "log")
  
  if (is.data.frame(counts_file) || is.matrix(counts_file)) {
    counts_data <- as.data.frame(counts_file)
  } else {
    # Determine format if auto
    if (counts_format == "auto") {
      ext <- tolower(tools::file_ext(counts_file))
      if (ext %in% c("csv")) {
        counts_format <- "csv"
      } else if (ext %in% c("tsv", "txt")) {
        counts_format <- "tsv"
      } else {
        stop("Could not automatically determine counts file format. Please specify explicitly.")
      }
    }
    
    # Read the file
    if (counts_format == "csv") {
      counts_data <- utils::read.csv(counts_file, row.names = 1, check.names = FALSE)
    } else if (counts_format %in% c("tsv", "txt")) {
      counts_data <- utils::read.delim(counts_file, row.names = 1, check.names = FALSE)
    } else {
      stop("Unsupported counts format: ", counts_format)
    }
  }
  
  # Import metadata
  tidyseq_obj <- add_message(tidyseq_obj, "Importing metadata", "log")
  
  if (is.data.frame(metadata_file)) {
    metadata_data <- metadata_file
  } else {
    # Determine format if auto
    if (metadata_format == "auto") {
      ext <- tolower(tools::file_ext(metadata_file))
      if (ext %in% c("csv")) {
        metadata_format <- "csv"
      } else if (ext %in% c("tsv", "txt")) {
        metadata_format <- "tsv"
      } else {
        stop("Could not automatically determine metadata file format. Please specify explicitly.")
      }
    }
    
    # Read the file
    if (metadata_format == "csv") {
      metadata_data <- utils::read.csv(metadata_file, check.names = FALSE)
    } else if (metadata_format %in% c("tsv", "txt")) {
      metadata_data <- utils::read.delim(metadata_file, check.names = FALSE)
    } else {
      stop("Unsupported metadata format: ", metadata_format)
    }
  }
  
  # Identify sample column in metadata if not specified
  if (is.null(sample_column)) {
    sample_column <- colnames(metadata_data)[1]
    tidyseq_obj <- add_message(
      tidyseq_obj, 
      paste0("Sample column not specified, using first column: '", sample_column, "'"),
      "log"
    )
  }
  
  # Convert metadata to data frame with rownames as sample IDs
  rownames(metadata_data) <- metadata_data[[sample_column]]
  
  # Validate matching of sample IDs
  if (validate_ids) {
    tidyseq_obj <- add_message(tidyseq_obj, "Validating sample IDs", "log")
    
    # Get sample IDs from counts data
    count_samples <- colnames(counts_data)
    
    # Get sample IDs from metadata
    metadata_samples <- rownames(metadata_data)
    
    # Check for missing samples in counts
    missing_in_counts <- setdiff(metadata_samples, count_samples)
    if (length(missing_in_counts) > 0) {
      warning_msg <- paste0(
        "The following samples are in metadata but missing in counts data: ",
        paste(missing_in_counts, collapse = ", ")
      )
      tidyseq_obj <- add_message(tidyseq_obj, warning_msg, "warning")
    }
    
    # Check for missing samples in metadata
    missing_in_metadata <- setdiff(count_samples, metadata_samples)
    if (length(missing_in_metadata) > 0) {
      warning_msg <- paste0(
        "The following samples are in counts data but missing in metadata: ",
        paste(missing_in_metadata, collapse = ", ")
      )
      tidyseq_obj <- add_message(tidyseq_obj, warning_msg, "warning")
    }
    
    # If no common samples, stop
    common_samples <- intersect(count_samples, metadata_samples)
    if (length(common_samples) == 0) {
      stop("No common samples found between counts data and metadata")
    }
    
    # Filter to common samples
    counts_data <- counts_data[, common_samples, drop = FALSE]
    metadata_data <- metadata_data[common_samples, , drop = FALSE]
  }
  
  # Store in tidyseq object
  tidyseq_obj$raw_counts <- counts_data
  tidyseq_obj$metadata <- metadata_data
  
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    paste0("Imported ", nrow(counts_data), " genes and ", ncol(counts_data), " samples"),
    "log"
  )
  
  return(tidyseq_obj)
}

#' Filter low-count genes from RNA-Seq data
#'
#' @description Remove genes with low counts across samples
#'
#' @param tidyseq_object A tidyseq object with raw counts
#' @param min_count Minimum count required in at least min_samples samples (default: 10)
#' @param min_samples Minimum number of samples that must have min_count (default: 2)
#' @param count_type Whether to filter on raw or CPM values (default: "raw")
#'
#' @return A tidyseq object with filtered counts
#'
#' @export
filter_low_counts <- function(tidyseq_object, 
                               min_count = 10, 
                               min_samples = 2,
                               count_type = "raw") {
  
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  if (is.null(tidyseq_object$raw_counts)) {
    stop("No raw counts found in tidyseq object. Run import_data() first.")
  }
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Filtering low-count genes (min_count=", min_count, 
           ", min_samples=", min_samples, ", count_type=", count_type, ")"),
    "log"
  )
  
  counts <- tidyseq_object$raw_counts
  
  # Convert to CPM if specified
  if (count_type == "cpm") {
    library(edgeR)
    cpm_counts <- edgeR::cpm(counts)
    
    # Filter based on CPM
    keep <- rowSums(cpm_counts >= min_count) >= min_samples
    
  } else if (count_type == "raw") {
    # Filter based on raw counts
    keep <- rowSums(counts >= min_count) >= min_samples
    
  } else {
    stop("Invalid count_type. Must be 'raw' or 'cpm'.")
  }
  
  # Apply filter
  filtered_counts <- counts[keep, , drop = FALSE]
  
  # Store results
  tidyseq_object$filtered_counts <- filtered_counts
  
  # Store parameters
  tidyseq_object$parameters$filtering <- list(
    min_count = min_count,
    min_samples = min_samples,
    count_type = count_type
  )
  
  # Add log message
  genes_removed <- nrow(counts) - nrow(filtered_counts)
  genes_kept <- nrow(filtered_counts)
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Removed ", genes_removed, " low-count genes (", 
           round(genes_removed/nrow(counts)*100, 1), "%). Kept ", 
           genes_kept, " genes."),
    "log"
  )
  
  return(tidyseq_object)
}

#' Normalize RNA-Seq count data
#'
#' @description Normalize RNA-Seq count data using various methods
#'
#' @param tidyseq_object A tidyseq object with raw or filtered counts
#' @param method Normalization method. One of "DESeq2", "TMM", "RLE", "CPM", or "TPM" (default: "DESeq2")
#' @param use_filtered Whether to use filtered counts (default: TRUE)
#' @param design Formula specifying the design for DESeq2 normalization (default: ~ 1)
#' @param gene_lengths Vector of gene lengths (required for TPM normalization)
#'
#' @return A tidyseq object with normalized counts
#'
#' @export
normalize_counts <- function(tidyseq_object,
                              method = "DESeq2",
                              use_filtered = TRUE,
                              design = ~ 1,
                              gene_lengths = NULL) {
  
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  # Check if counts are available
  if (use_filtered) {
    if (is.null(tidyseq_object$filtered_counts)) {
      stop("No filtered counts found in tidyseq object. Run filter_low_counts() first.")
    }
    counts <- tidyseq_object$filtered_counts
  } else {
    if (is.null(tidyseq_object$raw_counts)) {
      stop("No raw counts found in tidyseq object. Run import_data() first.")
    }
    counts <- tidyseq_object$raw_counts
  }
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Normalizing counts with method: ", method),
    "log"
  )
  
  # Normalize based on method
  if (method == "DESeq2") {
    # Check if metadata is available
    if (is.null(tidyseq_object$metadata)) {
      stop("Metadata is required for DESeq2 normalization. Run import_data() first.")
    }
    
    # Create DESeq2 object
    library(DESeq2)
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts,
      colData = tidyseq_object$metadata,
      design = design
    )
    
    # Estimate size factors
    dds <- DESeq2::estimateSizeFactors(dds)
    
    # Get normalized counts
    normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
    
    # Store DESeq2 object for later use
    tidyseq_object$deseq2_object <- dds
    
  } else if (method == "TMM") {
    # TMM normalization with edgeR
    library(edgeR)
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    
    # Get normalized counts
    normalized_counts <- edgeR::cpm(dge, log = FALSE) * 1e6
    
    # Store edgeR object for later use
    tidyseq_object$edger_object <- dge
    
  } else if (method == "RLE") {
    # RLE normalization with edgeR
    library(edgeR)
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "RLE")
    
    # Get normalized counts
    normalized_counts <- edgeR::cpm(dge, log = FALSE) * 1e6
    
  } else if (method == "CPM") {
    # CPM normalization with edgeR
    library(edgeR)
    normalized_counts <- edgeR::cpm(counts) * 1e6
    
  } else if (method == "TPM") {
    # Check if gene lengths are available
    if (is.null(gene_lengths)) {
      stop("Gene lengths are required for TPM normalization")
    }
    
    # Ensure gene lengths match count data
    if (!all(names(gene_lengths) %in% rownames(counts))) {
      stop("Gene names in gene_lengths do not match count data")
    }
    
    # Subset gene lengths to match count data
    gene_lengths <- gene_lengths[rownames(counts)]
    
    # Calculate TPM
    rate <- counts / gene_lengths
    normalized_counts <- t(t(rate) / colSums(rate) * 1e6)
    
  } else {
    stop("Unsupported normalization method: ", method)
  }
  
  # Store results
  tidyseq_object$normalized_counts <- normalized_counts
  tidyseq_object$normalization_method <- method
  
  # Store parameters
  tidyseq_object$parameters$normalization <- list(
    method = method,
    use_filtered = use_filtered
  )
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Normalized ", nrow(normalized_counts), " genes using ", method, " method"),
    "log"
  )
  
  return(tidyseq_object)
}