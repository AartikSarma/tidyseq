#' Run Gene Ontology (GO) enrichment analysis
#'
#' @description Perform Gene Ontology enrichment analysis on differential expression results
#'
#' @param tidyseq_object A tidyseq object with DE results
#' @param contrasts Contrasts to analyze (default: all)
#' @param ont GO ontology to test. One of "BP", "MF", "CC", or "ALL" (default: "ALL")
#' @param organism Organism for annotation. KEGG organism code or common name (default: "hsa")
#' @param gene_id_type Type of gene IDs. One of "ENSEMBL", "ENTREZID", "SYMBOL", etc. (default: "ENSEMBL")
#' @param p_cutoff P-value cutoff for DE genes (default: 0.05)
#' @param p_adj_cutoff Adjusted p-value cutoff for enrichment results (default: 0.05)
#' @param min_gs_size Minimum gene set size (default: 10)
#' @param max_gs_size Maximum gene set size (default: 500)
#'
#' @return A tidyseq object with GO enrichment results
#'
#' @export
run_go_enrichment <- function(tidyseq_object,
                               contrasts = NULL,
                               ont = "ALL",
                               organism = "hsa",
                               gene_id_type = "ENSEMBL",
                               p_cutoff = 0.05,
                               p_adj_cutoff = 0.05,
                               min_gs_size = 10,
                               max_gs_size = 500) {
  
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  # Check if DE results are available
  if (length(tidyseq_object$de_results) == 0) {
    stop("No differential expression results found. Run run_de_analysis() first.")
  }
  
  # Use all contrasts if not specified
  if (is.null(contrasts)) {
    contrasts <- names(tidyseq_object$de_results)
  } else {
    # Check if specified contrasts exist
    missing_contrasts <- setdiff(contrasts, names(tidyseq_object$de_results))
    if (length(missing_contrasts) > 0) {
      stop("Specified contrasts not found: ", paste(missing_contrasts, collapse = ", "))
    }
  }
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Running GO enrichment analysis for ", length(contrasts), " contrasts"),
    "log"
  )
  
  # Load required packages
  library(clusterProfiler)
  
  # Initialize results list if not present
  if (is.null(tidyseq_object$enrichment_results)) {
    tidyseq_object$enrichment_results <- list()
  }
  
  # Process each contrast
  for (contrast in contrasts) {
    tidyseq_object <- add_message(
      tidyseq_object, 
      paste0("Running GO enrichment for contrast: ", contrast),
      "log"
    )
    
    # Get DE results
    de_result <- tidyseq_object$de_results[[contrast]]
    
    # Get significant genes
    sig_genes <- de_result$gene_id[de_result$padj < p_cutoff]
    
    if (length(sig_genes) == 0) {
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("No significant genes found for contrast ", contrast, ". Skipping GO enrichment."),
        "warning"
      )
      next
    }
    
    # Get all genes as background
    all_genes <- de_result$gene_id
    
    # Convert gene IDs if needed
    if (gene_id_type != "ENSEMBL" && gene_id_type != "ENTREZID" && gene_id_type != "SYMBOL") {
      # Load appropriate annotation package
      if (organism == "hsa" || organism == "human") {
        library(org.Hs.eg.db)
        org_db <- org.Hs.eg.db
      } else if (organism == "mmu" || organism == "mouse") {
        library(org.Mm.eg.db)
        org_db <- org.Mm.eg.db
      } else if (organism == "rno" || organism == "rat") {
        library(org.Rn.eg.db)
        org_db <- org.Rn.eg.db
      } else {
        stop("Unsupported organism: ", organism)
      }
      
      # Convert gene IDs
      sig_genes <- clusterProfiler::bitr(sig_genes, fromType = gene_id_type, 
                                        toType = "ENTREZID", OrgDb = org_db)$ENTREZID
      all_genes <- clusterProfiler::bitr(all_genes, fromType = gene_id_type, 
                                        toType = "ENTREZID", OrgDb = org_db)$ENTREZID
    }
    
    # Run GO enrichment analysis
    go_results <- list()
    
    # Process each ontology
    ontologies <- if (ont == "ALL") c("BP", "MF", "CC") else ont
    
    for (current_ont in ontologies) {
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("Running ", current_ont, " ontology enrichment for contrast: ", contrast),
        "log"
      )
      
      # Run enrichment
      go_result <- clusterProfiler::enrichGO(
        gene = sig_genes,
        universe = all_genes,
        OrgDb = org_db,
        ont = current_ont,
        pAdjustMethod = "BH",
        pvalueCutoff = p_cutoff,
        qvalueCutoff = p_adj_cutoff,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      )
      
      # Store results
      go_results[[current_ont]] <- go_result
      
      # Log result summary
      if (length(go_result@result$ID) > 0) {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("Found ", length(go_result@result$ID), " enriched ", current_ont, " terms for contrast ", contrast),
          "log"
        )
      } else {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("No enriched ", current_ont, " terms found for contrast ", contrast),
          "log"
        )
      }
    }
    
    # Store all GO results for this contrast
    tidyseq_object$enrichment_results[[paste0(contrast, "_go")]] <- go_results
  }
  
  # Store enrichment parameters
  tidyseq_object$parameters$go_enrichment <- list(
    ont = ont,
    organism = organism,
    gene_id_type = gene_id_type,
    p_cutoff = p_cutoff,
    p_adj_cutoff = p_adj_cutoff,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )
  
  return(tidyseq_object)
}

#' Run KEGG pathway enrichment analysis
#'
#' @description Perform KEGG pathway enrichment analysis on differential expression results
#'
#' @param tidyseq_object A tidyseq object with DE results
#' @param contrasts Contrasts to analyze (default: all)
#' @param organism Organism for KEGG. KEGG organism code (default: "hsa")
#' @param gene_id_type Type of gene IDs. One of "ENSEMBL", "ENTREZID", "SYMBOL", etc. (default: "ENSEMBL")
#' @param p_cutoff P-value cutoff for DE genes (default: 0.05)
#' @param p_adj_cutoff Adjusted p-value cutoff for enrichment results (default: 0.05)
#' @param min_gs_size Minimum gene set size (default: 10)
#' @param max_gs_size Maximum gene set size (default: 500)
#'
#' @return A tidyseq object with KEGG enrichment results
#'
#' @export
run_kegg_enrichment <- function(tidyseq_object,
                                 contrasts = NULL,
                                 organism = "hsa",
                                 gene_id_type = "ENSEMBL",
                                 p_cutoff = 0.05,
                                 p_adj_cutoff = 0.05,
                                 min_gs_size = 10,
                                 max_gs_size = 500) {
  
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  # Check if DE results are available
  if (length(tidyseq_object$de_results) == 0) {
    stop("No differential expression results found. Run run_de_analysis() first.")
  }
  
  # Use all contrasts if not specified
  if (is.null(contrasts)) {
    contrasts <- names(tidyseq_object$de_results)
  } else {
    # Check if specified contrasts exist
    missing_contrasts <- setdiff(contrasts, names(tidyseq_object$de_results))
    if (length(missing_contrasts) > 0) {
      stop("Specified contrasts not found: ", paste(missing_contrasts, collapse = ", "))
    }
  }
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Running KEGG enrichment analysis for ", length(contrasts), " contrasts"),
    "log"
  )
  
  # Load required packages
  library(clusterProfiler)
  
  # Initialize results list if not present
  if (is.null(tidyseq_object$enrichment_results)) {
    tidyseq_object$enrichment_results <- list()
  }
  
  # Process each contrast
  for (contrast in contrasts) {
    tidyseq_object <- add_message(
      tidyseq_object, 
      paste0("Running KEGG enrichment for contrast: ", contrast),
      "log"
    )
    
    # Get DE results
    de_result <- tidyseq_object$de_results[[contrast]]
    
    # Get significant genes
    sig_genes <- de_result$gene_id[de_result$padj < p_cutoff]
    
    if (length(sig_genes) == 0) {
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("No significant genes found for contrast ", contrast, ". Skipping KEGG enrichment."),
        "warning"
      )
      next
    }
    
    # Get all genes as background
    all_genes <- de_result$gene_id
    
    # Convert gene IDs if needed
    if (gene_id_type != "ENTREZID") {
      # Load appropriate annotation package
      if (organism == "hsa" || organism == "human") {
        library(org.Hs.eg.db)
        org_db <- org.Hs.eg.db
      } else if (organism == "mmu" || organism == "mouse") {
        library(org.Mm.eg.db)
        org_db <- org.Mm.eg.db
      } else if (organism == "rno" || organism == "rat") {
        library(org.Rn.eg.db)
        org_db <- org.Rn.eg.db
      } else {
        stop("Unsupported organism: ", organism)
      }
      
      # Convert gene IDs
      sig_genes <- clusterProfiler::bitr(sig_genes, fromType = gene_id_type, 
                                        toType = "ENTREZID", OrgDb = org_db)$ENTREZID
      all_genes <- clusterProfiler::bitr(all_genes, fromType = gene_id_type, 
                                        toType = "ENTREZID", OrgDb = org_db)$ENTREZID
    }
    
    # Run KEGG enrichment analysis
    kegg_result <- clusterProfiler::enrichKEGG(
      gene = sig_genes,
      universe = all_genes,
      organism = organism,
      pvalueCutoff = p_cutoff,
      pAdjustMethod = "BH",
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      keyType = "ncbi-geneid"
    )
    
    # Store results
    tidyseq_object$enrichment_results[[paste0(contrast, "_kegg")]] <- kegg_result
    
    # Log result summary
    if (length(kegg_result@result$ID) > 0) {
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("Found ", length(kegg_result@result$ID), " enriched KEGG pathways for contrast ", contrast),
        "log"
      )
    } else {
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("No enriched KEGG pathways found for contrast ", contrast),
        "log"
      )
    }
  }
  
  # Store enrichment parameters
  tidyseq_object$parameters$kegg_enrichment <- list(
    organism = organism,
    gene_id_type = gene_id_type,
    p_cutoff = p_cutoff,
    p_adj_cutoff = p_adj_cutoff,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )
  
  return(tidyseq_object)
}

#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' @description Perform Gene Set Enrichment Analysis on ranked gene lists
#'
#' @param tidyseq_object A tidyseq object with DE results
#' @param contrasts Contrasts to analyze (default: all)
#' @param gene_sets Gene sets to test. One of "GO", "KEGG", or a custom GMT file (default: "GO")
#' @param organism Organism for annotation (default: "hsa")
#' @param gene_id_type Type of gene IDs (default: "ENSEMBL")
#' @param rank_by Method to rank genes. One of "pvalue", "log2fc", or "stat" (default: "stat")
#' @param p_adj_cutoff Adjusted p-value cutoff (default: 0.05)
#' @param min_gs_size Minimum gene set size (default: 10)
#' @param max_gs_size Maximum gene set size (default: 500)
#'
#' @return A tidyseq object with GSEA results
#'
#' @export
run_gsea <- function(tidyseq_object,
                      contrasts = NULL,
                      gene_sets = "GO",
                      organism = "hsa",
                      gene_id_type = "ENSEMBL",
                      rank_by = "stat",
                      p_adj_cutoff = 0.05,
                      min_gs_size = 10,
                      max_gs_size = 500) {
  
  if (!inherits(tidyseq_object, "tidyseq")) {
    stop("Object must be of class 'tidyseq'")
  }
  
  # Check if DE results are available
  if (length(tidyseq_object$de_results) == 0) {
    stop("No differential expression results found. Run run_de_analysis() first.")
  }
  
  # Use all contrasts if not specified
  if (is.null(contrasts)) {
    contrasts <- names(tidyseq_object$de_results)
  } else {
    # Check if specified contrasts exist
    missing_contrasts <- setdiff(contrasts, names(tidyseq_object$de_results))
    if (length(missing_contrasts) > 0) {
      stop("Specified contrasts not found: ", paste(missing_contrasts, collapse = ", "))
    }
  }
  
  tidyseq_object <- add_message(
    tidyseq_object, 
    paste0("Running GSEA for ", length(contrasts), " contrasts"),
    "log"
  )
  
  # Load required packages
  library(clusterProfiler)
  
  # Initialize results list if not present
  if (is.null(tidyseq_object$enrichment_results)) {
    tidyseq_object$enrichment_results <- list()
  }
  
  # Process each contrast
  for (contrast in contrasts) {
    tidyseq_object <- add_message(
      tidyseq_object, 
      paste0("Running GSEA for contrast: ", contrast),
      "log"
    )
    
    # Get DE results
    de_result <- tidyseq_object$de_results[[contrast]]
    
    # Create ranked gene list
    if (rank_by == "pvalue") {
      # Rank by -log10(pvalue) * sign(log2FoldChange)
      gene_list <- -log10(de_result$pvalue) * sign(de_result$log2FoldChange)
    } else if (rank_by == "log2fc") {
      # Rank by log2FoldChange
      gene_list <- de_result$log2FoldChange
    } else if (rank_by == "stat") {
      # Rank by stat
      if ("stat" %in% colnames(de_result)) {
        gene_list <- de_result$stat
      } else {
        # Calculate stat as log2FoldChange / lfcSE
        gene_list <- de_result$log2FoldChange / de_result$lfcSE
      }
    } else {
      stop("Invalid rank_by parameter. Must be one of: pvalue, log2fc, stat")
    }
    
    # Name the gene list
    names(gene_list) <- de_result$gene_id
    
    # Sort the list in decreasing order
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    # Remove any NA values
    gene_list <- gene_list[!is.na(gene_list)]
    
    # Convert gene IDs if needed
    if (gene_id_type != "ENTREZID" && gene_sets %in% c("GO", "KEGG")) {
      # Load appropriate annotation package
      if (organism == "hsa" || organism == "human") {
        library(org.Hs.eg.db)
        org_db <- org.Hs.eg.db
      } else if (organism == "mmu" || organism == "mouse") {
        library(org.Mm.eg.db)
        org_db <- org.Mm.eg.db
      } else if (organism == "rno" || organism == "rat") {
        library(org.Rn.eg.db)
        org_db <- org.Rn.eg.db
      } else {
        stop("Unsupported organism: ", organism)
      }
      
      # Convert gene IDs
      gene_ids <- names(gene_list)
      converted_ids <- clusterProfiler::bitr(gene_ids, fromType = gene_id_type, 
                                           toType = "ENTREZID", OrgDb = org_db)
      
      # Match original list with converted IDs
      gene_list <- gene_list[gene_ids %in% converted_ids[[gene_id_type]]]
      names(gene_list) <- converted_ids$ENTREZID[match(names(gene_list), converted_ids[[gene_id_type]])]
    }
    
    # Run GSEA based on gene sets
    if (gene_sets == "GO") {
      # Run GO GSEA
      go_results <- list()
      
      # Process each ontology
      for (ont in c("BP", "MF", "CC")) {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("Running ", ont, " GSEA for contrast: ", contrast),
          "log"
        )
        
        # Run enrichment
        gsea_result <- clusterProfiler::gseGO(
          geneList = gene_list,
          OrgDb = org_db,
          ont = ont,
          minGSSize = min_gs_size,
          maxGSSize = max_gs_size,
          pvalueCutoff = p_adj_cutoff,
          pAdjustMethod = "BH",
          verbose = FALSE
        )
        
        # Store results
        go_results[[ont]] <- gsea_result
        
        # Log result summary
        if (nrow(gsea_result@result) > 0) {
          tidyseq_object <- add_message(
            tidyseq_object, 
            paste0("Found ", nrow(gsea_result@result), " enriched ", ont, " gene sets for contrast ", contrast),
            "log"
          )
        } else {
          tidyseq_object <- add_message(
            tidyseq_object, 
            paste0("No enriched ", ont, " gene sets found for contrast ", contrast),
            "log"
          )
        }
      }
      
      # Store all GO GSEA results for this contrast
      tidyseq_object$enrichment_results[[paste0(contrast, "_gsea_go")]] <- go_results
      
    } else if (gene_sets == "KEGG") {
      # Run KEGG GSEA
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("Running KEGG GSEA for contrast: ", contrast),
        "log"
      )
      
      gsea_result <- clusterProfiler::gseKEGG(
        geneList = gene_list,
        organism = organism,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        pvalueCutoff = p_adj_cutoff,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
      
      # Store results
      tidyseq_object$enrichment_results[[paste0(contrast, "_gsea_kegg")]] <- gsea_result
      
      # Log result summary
      if (nrow(gsea_result@result) > 0) {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("Found ", nrow(gsea_result@result), " enriched KEGG pathways for contrast ", contrast),
          "log"
        )
      } else {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("No enriched KEGG pathways found for contrast ", contrast),
          "log"
        )
      }
      
    } else {
      # Custom gene sets from GMT file
      tidyseq_object <- add_message(
        tidyseq_object, 
        paste0("Running custom GSEA for contrast: ", contrast, " using ", gene_sets),
        "log"
      )
      
      # Check if file exists
      if (!file.exists(gene_sets)) {
        stop("Custom gene sets file not found: ", gene_sets)
      }
      
      # Read GMT file
      custom_gs <- clusterProfiler::read.gmt(gene_sets)
      
      # Run GSEA
      gsea_result <- clusterProfiler::GSEA(
        geneList = gene_list,
        TERM2GENE = custom_gs,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        pvalueCutoff = p_adj_cutoff,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
      
      # Store results
      tidyseq_object$enrichment_results[[paste0(contrast, "_gsea_custom")]] <- gsea_result
      
      # Log result summary
      if (nrow(gsea_result@result) > 0) {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("Found ", nrow(gsea_result@result), " enriched gene sets for contrast ", contrast),
          "log"
        )
      } else {
        tidyseq_object <- add_message(
          tidyseq_object, 
          paste0("No enriched gene sets found for contrast ", contrast),
          "log"
        )
      }
    }
  }
  
  # Store GSEA parameters
  tidyseq_object$parameters$gsea <- list(
    gene_sets = gene_sets,
    organism = organism,
    gene_id_type = gene_id_type,
    rank_by = rank_by,
    p_adj_cutoff = p_adj_cutoff,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size
  )
  
  return(tidyseq_object)
}