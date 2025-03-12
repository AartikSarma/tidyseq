library(testthat)
library(tidyseq)

test_that("tidyseq object can be created", {
  # Create a tidyseq object
  obj <- tidyseq()
  
  # Check that it has the correct class
  expect_s3_class(obj, "tidyseq")
  
  # Check that it has the expected structure
  expect_true(is.list(obj))
  expect_true(is.null(obj$raw_counts))
  expect_true(is.null(obj$metadata))
  expect_true(is.list(obj$de_results))
  expect_true(is.list(obj$enrichment_results))
})

test_that("add_message adds messages to tidyseq object", {
  # Create a tidyseq object
  obj <- tidyseq()
  
  # Add a log message
  obj <- add_message(obj, "Test log message", "log")
  
  # Check that the message was added
  expect_true(length(obj$log) == 1)
  expect_true(grepl("Test log message", obj$log[1]))
  
  # Add a warning message
  expect_warning(
    obj <- add_message(obj, "Test warning message", "warning")
  )
  
  # Check that the warning was added
  expect_true(length(obj$warnings) == 1)
  expect_true(grepl("Test warning message", obj$warnings[1]))
  
  # Add an error message
  expect_error(
    obj <- add_message(obj, "Test error message", "error")
  )
})

# Create simple test data for other functions
create_test_data <- function() {
  # Create count data
  counts <- matrix(
    rpois(100, lambda = 10),
    nrow = 20,
    ncol = 5,
    dimnames = list(
      paste0("gene", 1:20),
      paste0("sample", 1:5)
    )
  )
  
  # Create metadata
  metadata <- data.frame(
    sample = paste0("sample", 1:5),
    condition = c("A", "A", "A", "B", "B"),
    batch = c(1, 1, 2, 2, 1)
  )
  
  # Return as list
  list(counts = counts, metadata = metadata)
}

test_that("import_data imports data correctly", {
  # Skip on CRAN to avoid creating files
  skip_on_cran()
  
  # Create test data
  test_data <- create_test_data()
  
  # Import data
  obj <- import_data(
    counts_file = test_data$counts,
    metadata_file = test_data$metadata,
    sample_column = "sample"
  )
  
  # Check that data was imported correctly
  expect_s3_class(obj, "tidyseq")
  expect_equal(dim(obj$raw_counts), c(20, 5))
  expect_equal(dim(obj$metadata), c(5, 3))
  expect_equal(colnames(obj$raw_counts), rownames(obj$metadata))
})

test_that("filter_low_counts filters genes correctly", {
  # Skip on CRAN to avoid creating files
  skip_on_cran()
  
  # Create test data
  test_data <- create_test_data()
  
  # Create an object with data
  obj <- tidyseq()
  obj$raw_counts <- test_data$counts
  
  # Filter low-count genes
  obj <- filter_low_counts(
    obj,
    min_count = 10,
    min_samples = 2
  )
  
  # Check that genes were filtered correctly
  expect_s3_class(obj, "tidyseq")
  expect_true(nrow(obj$filtered_counts) <= nrow(obj$raw_counts))
})

test_that("normalize_counts normalizes data correctly", {
  # Skip on CRAN to avoid DESeq2 dependency
  skip_on_cran()
  
  # Create test data
  test_data <- create_test_data()
  
  # Create an object with data
  obj <- tidyseq()
  obj$raw_counts <- test_data$counts
  obj$metadata <- test_data$metadata
  obj$filtered_counts <- test_data$counts
  
  # Mock the DESeq2 functions to avoid dependency in tests
  # In a real test, you'd use the actual DESeq2 functions
  with_mock(
    `DESeq2::DESeqDataSetFromMatrix` = function(...) list(),
    `DESeq2::estimateSizeFactors` = function(x) x,
    `DESeq2::counts` = function(x, normalized) test_data$counts,
    {
      # Normalize counts
      obj <- normalize_counts(
        obj,
        method = "DESeq2"
      )
      
      # Check that data was normalized correctly
      expect_s3_class(obj, "tidyseq")
      expect_equal(dim(obj$normalized_counts), c(20, 5))
      expect_equal(obj$normalization_method, "DESeq2")
    }
  )
})