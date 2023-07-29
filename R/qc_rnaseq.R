#' Quality control and exploratory analysis for RNAseq count data
#'
#' This function performs basic quality control and exploratory analysis 
#' for RNAseq count data. It takes as input a metadata tibble and a counts tibble
#' and confirms that the sample IDs match between the two. It renames the gene ID
#' column in the counts tibble, checks that the remaining columns are numeric, 
#' and calculates QC metrics including total counts, nonzero counts, and 
#' ribosomal/mitochondrial percentages per sample. It generates plots for these
#' QC metrics using ggplot2. It also generates a PCA plot and cluster dendrogram 
#' for the samples.
#'
#' @param metadata A tibble with sample metadata
#' @param counts A tibble with RNAseq counts, with samples as columns
#' @param sample_col Name of sample ID column in metadata (default 'sampleid')
#' @param gene_col Name of gene ID column in counts (default 'ensembl_gene_id')
#'
#' @return A list with two elements:
#' \itemize{
#'   \item{qc_metrics}{A tibble with QC metrics per sample}
#'   \item{plots}{A list containing QC metric plots and PCA and clustering plots}
#' }
#'
#' @examples
#' metadata <- tibble(sampleid = c("s1", "s2", "s3"), 
#'                    tissue = c("liver", "brain", "kidney"))
#' counts <- tibble(ensembl_gene_id = c("ENSG1", "ENSG2"),
#'                  s1 = c(10, 20),
#'                  s2 = c(15, 50), 
#'                  s3 = c(20, 100))
#' qc_data <- rnaseq_qc(metadata, counts)
#'
#' @export
rnaseq_qc <- function(metadata, counts, 
                      sample_col = "sampleid",
                      gene_col = "ensembl_gene_id") {
  
  # Check that sample IDs match between metadata and counts
  meta_samples <- unique(metadata[[sample_col]])
  count_samples <- colnames(counts)[-1]
  if (!all(meta_samples %in% count_samples)) {
    stop("Sample IDs in metadata do not match count columns") 
  }
  
  # Check that counts data is appropriate
  if(!is.character(counts[[1]])) {
    stop("Counts data does not contain a character column for gene IDs")
  }
  if(sum(sapply(counts[,2:ncol(counts)], is.numeric)) != ncol(counts)-1) {
    stop("Counts data contains non-numeric columns")
  }
  
  # Rename gene ID column
  colnames(counts)[1] <- "gene_id" 
  
  # Calculate QC metrics
  counts <- mutate(counts, total = rowSums(.[,-1]))
  qc_metrics <- counts %>%
    pivot_longer(cols = -gene_id, names_to = "sample", values_to = "count") %>%
    left_join(metadata, by = c("sample" = sample_col)) %>%
    group_by(sample) %>%
    summarise(total = sum(count), 
              zeros = sum(count == 0),
              pct_ribo = mean(grepl("RPS|RPL", gene_id)) * 100, 
              pct_mito = mean(grepl("^MT-", gene_id)) * 100)
  
  # PCA plot
  vst_counts <- counts %>% 
    column_to_rownames("gene_id") %>%
    t() %>%
    as.matrix() %>%
    vst() %>%
    t() %>%
    as_tibble()
  pca_data <- prcomp(t(scale(t(vst_counts[,2:ncol(vst_counts)]))), 
                     center = TRUE)
  pca_plot <- ggplot(as_tibble(pca_data$x)) +
    aes(PC1, PC2, color = metadata[[sample_col]]) +
    geom_point() +
    labs(title="PCA",
         x = paste0("PC1 (",scales::percent(pca_data$sdev[1]^2/sum(pca_data$sdev^2)), ")"),
         y = paste0("PC2 (",scales::percent(pca_data$sdev[2]^2/sum(pca_data$sdev^2)), ")"))
  
  # Clustering
  dist_mat <- dist(t(vst_counts[,2:ncol(vst_counts)]), method = "euclidean")
  hclust_data <- hclust(dist_mat)
  clust_plot <- ggplot(as_tibble(hclust_data$height)) +
    aes(x = hclust_data$height) +
    geom_histogram(bins=30) +
    labs(x = "Height", y = "Count")
  
  # Return results
  list(qc_metrics = qc_metrics,
       plots = list(qc = qc_plots,
                    pca = pca_plot,
                    clustering = clust_plot))
  
}