# TidySeq: User-Friendly RNA-Seq Analysis for Researchers with Limited Programming Experience

TidySeq is an R package designed to simplify RNA-Seq analysis for researchers with limited programming experience. It provides a streamlined workflow for data import, quality control, normalization, differential expression analysis, pathway enrichment, and visualization with robust error handling and guidance.

## Features

- **User-Friendly Interface**: Simple functions with sensible defaults for common RNA-Seq analyses
- **Comprehensive Workflow**: Complete pipeline from raw counts to functional enrichment
- **Informative Visualization**: Publication-ready plots for data exploration and result interpretation
- **Robust Error Handling**: Clear error messages and validation to prevent common mistakes
- **Automated Reporting**: Generate comprehensive HTML reports with a single function
- **Educational Explanations**: Built-in explanations with example R code to help you understand what each function is doing
- **Flexible Design**: Use the complete workflow or individual components as needed

## Installation

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install TidySeq from GitHub
devtools::install_github("yourusername/tidyseq")
```

## Quick Start

```r
library(tidyseq)

# Run complete analysis workflow
results <- run_rnaseq_workflow(
  counts_file = "counts.csv",
  metadata_file = "metadata.csv",
  condition_column = "treatment",
  output_dir = "results"
)

# Print results summary
print(results)
```

## Step-by-Step Analysis

```r
# Create a tidyseq object
rnaseq <- tidyseq()

# Import data
rnaseq <- import_data(
  rnaseq,
  counts_file = "counts.csv",
  metadata_file = "metadata.csv"
)

# Filter low-count genes
rnaseq <- filter_low_counts(
  rnaseq,
  min_count = 10,
  min_samples = 2
)

# Normalize counts
rnaseq <- normalize_counts(
  rnaseq,
  method = "DESeq2"
)

# Create PCA plot
pca_plot <- plot_pca(
  rnaseq,
  color_by = "treatment"
)

# Run differential expression analysis
rnaseq <- create_de_contrasts(
  rnaseq,
  condition_column = "treatment"
)

rnaseq <- run_de_analysis(
  rnaseq,
  method = "DESeq2"
)

# Create volcano plot
volcano_plot <- plot_volcano(
  rnaseq,
  contrast = "treatment_treated_vs_control"
)

# Run enrichment analysis
rnaseq <- run_go_enrichment(
  rnaseq,
  organism = "hsa"
)

# Generate comprehensive report
rnaseq <- create_full_report(
  rnaseq,
  output_dir = "results"
)
```

## Key Components

1. **Data Import & Preparation**
   - Import count data and metadata
   - Filter low-count genes
   - Normalize count data

2. **Quality Control**
   - Sample correlation analysis
   - Principal component analysis
   - Count distribution visualization

3. **Differential Expression Analysis**
   - Multiple DE methods (DESeq2, limma-voom, edgeR)
   - Flexible contrast specification
   - Interactive volcano and MA plots

4. **Functional Enrichment**
   - Gene Ontology (GO) enrichment
   - KEGG pathway analysis
   - Gene Set Enrichment Analysis (GSEA)

5. **Visualization & Reporting**
   - Publication-ready plots
   - Interactive visualizations
   - Comprehensive HTML reports

## Educational Explanations

TidySeq includes a unique feature that explains what each function is doing with example R code, helping you learn how RNA-Seq analysis works:

```r
# By default, explanations are enabled
tidyseq_explanations_enabled()  # Returns TRUE

# Run an analysis with explanations
rnaseq <- import_data(
  counts_file = "counts.csv",
  metadata_file = "metadata.csv"
)
# This will display an explanation of what import_data() does

# Turn explanations off globally
tidyseq_explanations_off()

# Turn explanations on globally
tidyseq_explanations_on()

# Control explanations for a specific function call
rnaseq <- filter_low_counts(
  rnaseq,
  min_count = 10,
  min_samples = 2,
  explain = TRUE  # Override global setting for this call
)
```

## Documentation

For more information, please refer to the package vignettes:

```r
browseVignettes("tidyseq")
```

## Contribution

Contributions to TidySeq are welcome! Please feel free to submit a pull request or open an issue on GitHub.

## License

This project is licensed under the MIT License - see the LICENSE file for details.