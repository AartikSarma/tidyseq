#' PCA plot for RNAseq count data 
#'
#' Generates a PCA plot colored by a sample annotation
#'
#' @param metadata A tibble with sample metadata 
#' @param counts A tibble with RNAseq count data
#' @param annotation Column in metadata to use for annotation (default NULL)
#'
#' @return A ggplot2 PCA plot
#' 
#' @examples
#' metadata <- tibble(sample = c("s1", "s2", "s3"),
#'                    tissue = c("liver", "brain", "kidney"))  
#' counts <- tibble(gene_id = c("g1", "g2"),
#'                  s1 = c(10, 20), 
#'                  s2 = c(15, 50),
#'                  s3 = c(20, 100))
#' pca_plot(metadata, counts, annotation = "tissue")
#'
#' @export

pca_plot <- function(metadata, counts, annotation = NULL) {
  
  # Check annotation column
  if(!is.null(annotation)) {
    if(!annotation %in% colnames(metadata)) {
      stop("Annotation column not found in metadata")
    }
  }
  
  # Transform counts
  vst_counts <- counts %>%
    column_to_rownames("gene_id") %>% 
    t() %>%
    as.matrix() %>%
    vst() %>%
    t() %>%
    as_tibble()
  
  # Run PCA
  pca_data <- prcomp(t(scale(t(vst_counts[,2:ncol(vst_counts)]))),
                     center = TRUE)
  
  # Generate plot
  plot <- ggplot(as_tibble(pca_data$x)) +
    aes(PC1, PC2) +
    geom_point(aes_string(color = annotation)) +
    labs(title = "PCA Plot",
         x = paste0("PC1 (",scales::percent(pca_data$sdev[1]^2/sum(pca_data$sdev^2)), ")"),
         y = paste0("PC2 (",scales::percent(pca_data$sdev[2]^2/sum(pca_data$sdev^2)), ")"))
  
  return(plot)
}