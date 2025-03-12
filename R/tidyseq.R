#' TidySeq: User-Friendly RNA-Seq Analysis for Researchers with Limited Programming Experience
#'
#' @description TidySeq simplifies RNA-Seq analysis for researchers with limited 
#' programming experience. It provides a streamlined workflow for data import, 
#' quality control, normalization, differential expression analysis, pathway 
#' enrichment, and visualization with robust error handling and guidance.
#'
#' @docType package
#' @name tidyseq-package
NULL

#' Complete RNA-Seq analysis workflow in one function
#'
#' @description Run the complete RNA-Seq analysis workflow from raw counts to enrichment analysis
#'
#' @param counts_file Path to count data file or a data frame/matrix of counts
#' @param metadata_file Path to metadata file or a data frame of metadata
#' @param condition_column Column in metadata containing condition information for DE analysis
#' @param reference_level Reference level for condition comparisons (optional)
#' @param output_dir Output directory for results and reports
#' @param filter_params List of parameters for filtering (passed to filter_low_counts)
#' @param norm_params List of parameters for normalization (passed to normalize_counts)
#' @param de_params List of parameters for DE analysis (passed to run_de_analysis)
#' @param enrich_params List of parameters for enrichment analysis (optional)
#' @param create_reports Whether to create HTML reports (default: TRUE)
#'
#' @return A tidyseq object with all analysis results
#'
#' @examples
#' \dontrun{
#' # Basic workflow with default parameters
#' results <- run_rnaseq_workflow(
#'   counts_file = "counts.csv",
#'   metadata_file = "metadata.csv",
#'   condition_column = "treatment",
#'   output_dir = "results"
#' )
#'
#' # Advanced workflow with custom parameters
#' results <- run_rnaseq_workflow(
#'   counts_file = "counts.csv",
#'   metadata_file = "metadata.csv",
#'   condition_column = "treatment",
#'   reference_level = "control",
#'   output_dir = "results",
#'   filter_params = list(min_count = 5, min_samples = 3, count_type = "cpm"),
#'   norm_params = list(method = "TMM"),
#'   de_params = list(method = "limma-voom", alpha = 0.1),
#'   enrich_params = list(run_go = TRUE, run_kegg = TRUE, run_gsea = FALSE)
#' )
#' }
#'
#' @export
run_rnaseq_workflow <- function(counts_file,
                                 metadata_file,
                                 condition_column,
                                 reference_level = NULL,
                                 output_dir = "tidyseq_results",
                                 filter_params = list(),
                                 norm_params = list(),
                                 de_params = list(),
                                 enrich_params = list(),
                                 create_reports = TRUE) {
  
  # Start TidySeq analysis
  tidyseq_obj <- tidyseq()
  
  # Add initial log message
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "Starting TidySeq analysis workflow",
    "log"
  )
  
  # Step 1: Import data
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "STEP 1: Importing count data and metadata",
    "log"
  )
  
  tidyseq_obj <- import_data(
    counts_file = counts_file,
    metadata_file = metadata_file
  )
  
  # Step 2: Filter low-count genes
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "STEP 2: Filtering low-count genes",
    "log"
  )
  
  # Set default filtering parameters
  filter_defaults <- list(
    min_count = 10, 
    min_samples = 2,
    count_type = "raw"
  )
  
  # Update with user-provided parameters
  filter_args <- modifyList(filter_defaults, filter_params)
  
  # Run filtering
  tidyseq_obj <- do.call(filter_low_counts, c(list(tidyseq_object = tidyseq_obj), filter_args))
  
  # Step 3: Normalize counts
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "STEP 3: Normalizing count data",
    "log"
  )
  
  # Set default normalization parameters
  norm_defaults <- list(
    method = "DESeq2", 
    use_filtered = TRUE
  )
  
  # Update with user-provided parameters
  norm_args <- modifyList(norm_defaults, norm_params)
  
  # Run normalization
  tidyseq_obj <- do.call(normalize_counts, c(list(tidyseq_object = tidyseq_obj), norm_args))
  
  # Step 4: Create contrasts for DE analysis
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "STEP 4: Creating DE contrasts",
    "log"
  )
  
  tidyseq_obj <- create_de_contrasts(
    tidyseq_object = tidyseq_obj,
    condition_column = condition_column,
    reference_level = reference_level
  )
  
  # Step 5: Run DE analysis
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "STEP 5: Running differential expression analysis",
    "log"
  )
  
  # Set default DE parameters
  de_defaults <- list(
    method = "DESeq2", 
    alpha = 0.05,
    lfc_threshold = 0,
    use_filtered = TRUE,
    use_normalized = TRUE
  )
  
  # Update with user-provided parameters
  de_args <- modifyList(de_defaults, de_params)
  
  # Run DE analysis
  tidyseq_obj <- do.call(run_de_analysis, c(list(tidyseq_object = tidyseq_obj), de_args))
  
  # Step 6: Run enrichment analysis
  # Check if enrichment analysis is requested
  run_enrichment <- is.null(enrich_params$run_enrichment) || enrich_params$run_enrichment
  run_go <- is.null(enrich_params$run_go) || enrich_params$run_go
  run_kegg <- is.null(enrich_params$run_kegg) || enrich_params$run_kegg
  run_gsea <- is.null(enrich_params$run_gsea) || enrich_params$run_gsea
  
  if (run_enrichment) {
    tidyseq_obj <- add_message(
      tidyseq_obj, 
      "STEP 6: Running functional enrichment analysis",
      "log"
    )
    
    # Set default enrichment parameters
    enrich_defaults <- list(
      organism = "hsa",
      gene_id_type = "ENSEMBL",
      p_cutoff = 0.05,
      p_adj_cutoff = 0.05
    )
    
    # Update with user-provided parameters
    enrich_args <- modifyList(enrich_defaults, enrich_params)
    
    # Run GO enrichment if requested
    if (run_go) {
      tidyseq_obj <- add_message(
        tidyseq_obj, 
        "Running GO enrichment analysis",
        "log"
      )
      
      tidyseq_obj <- do.call(run_go_enrichment, c(list(tidyseq_object = tidyseq_obj), enrich_args))
    }
    
    # Run KEGG enrichment if requested
    if (run_kegg) {
      tidyseq_obj <- add_message(
        tidyseq_obj, 
        "Running KEGG pathway enrichment analysis",
        "log"
      )
      
      tidyseq_obj <- do.call(run_kegg_enrichment, c(list(tidyseq_object = tidyseq_obj), enrich_args))
    }
    
    # Run GSEA if requested
    if (run_gsea) {
      tidyseq_obj <- add_message(
        tidyseq_obj, 
        "Running Gene Set Enrichment Analysis (GSEA)",
        "log"
      )
      
      tidyseq_obj <- do.call(run_gsea, c(list(tidyseq_object = tidyseq_obj), enrich_args))
    }
  }
  
  # Step 7: Create reports
  if (create_reports) {
    tidyseq_obj <- add_message(
      tidyseq_obj, 
      "STEP 7: Creating analysis reports",
      "log"
    )
    
    # Create directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create QC report
    tidyseq_obj <- create_qc_report(
      tidyseq_object = tidyseq_obj,
      output_dir = output_dir
    )
    
    # Create DE report
    tidyseq_obj <- create_de_report(
      tidyseq_object = tidyseq_obj,
      output_dir = output_dir
    )
    
    # Create enrichment report if enrichment was performed
    if (run_enrichment && length(tidyseq_obj$enrichment_results) > 0) {
      tidyseq_obj <- create_enrichment_report(
        tidyseq_object = tidyseq_obj,
        output_dir = output_dir
      )
    }
    
    # Create full report
    tidyseq_obj <- create_full_report(
      tidyseq_object = tidyseq_obj,
      output_dir = output_dir
    )
  }
  
  # Add final log message
  tidyseq_obj <- add_message(
    tidyseq_obj, 
    "TidySeq analysis workflow completed successfully",
    "log"
  )
  
  return(tidyseq_obj)
}

#' Save TidySeq analysis results
#'
#' @description Save a TidySeq object to a file for later use
#'
#' @param tidyseq_object A tidyseq object
#' @param file_path Path to save the object
#'
#' @return Invisibly returns the file path
#'
#' @export
save_analysis <- function(tidyseq_object, file_path) {
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  # Create directory if it doesn't exist
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Save the object
  saveRDS(tidyseq_object, file = file_path)
  
  # Log message
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("TidySeq analysis saved to: ", file_path),
    "log"
  )
  
  return(invisible(file_path))
}

#' Load TidySeq analysis results
#'
#' @description Load a TidySeq object from a file
#'
#' @param file_path Path to the saved tidyseq object
#'
#' @return A tidyseq object
#'
#' @export
load_analysis <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Load the object
  tidyseq_object <- readRDS(file_path)
  
  # Verify it's a tidyseq object
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Loaded object is not of class 'tidyseq'")
  }
  
  # Add log message
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("TidySeq analysis loaded from: ", file_path),
    "log"
  )
  
  return(tidyseq_object)
}