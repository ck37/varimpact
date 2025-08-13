# Test single column support in varimpact
# This addresses issue #6: "Single column in data"

library(testthat)
library(varimpact)

context("Single column support")

test_that("separate_factors_numerics handles single columns", {
  # Test single numeric column
  data_numeric <- data.frame(x1 = c(1.1, 2.2, 3.3, 4.4, 5.5))
  result <- separate_factors_numerics(data_numeric)
  expect_equal(ncol(result$df_factors), 0)
  expect_equal(ncol(result$df_numerics), 1)
  expect_equal(colnames(result$df_numerics), "x1")
  
  # Test single factor column
  data_factor <- data.frame(x1 = factor(c("A", "B", "C", "A", "B")))
  result <- separate_factors_numerics(data_factor)
  expect_equal(ncol(result$df_factors), 1)
  expect_equal(ncol(result$df_numerics), 0)
  expect_equal(colnames(result$df_factors), "x1")
  
  # Test single character column (should be converted to factor)
  data_char <- data.frame(x1 = c("A", "B", "C", "A", "B"), stringsAsFactors = FALSE)
  result <- separate_factors_numerics(data_char)
  expect_equal(ncol(result$df_factors), 1)
  expect_equal(ncol(result$df_numerics), 0)
  expect_equal(colnames(result$df_factors), "x1")
  expect_true(is.factor(result$df_factors$x1))
})

test_that("varimpact handles vector input", {
  # Create test data
  set.seed(1)
  N <- 50
  X_vector <- rnorm(N)
  Y <- rbinom(N, 1, plogis(0.2 * X_vector))
  
  # Test that vector input is converted to data frame and processed
  expect_warning(
    vim <- varimpact(Y = Y, data = X_vector, verbose = FALSE, V = 2L),
    "Using single variable for variable importance analysis"
  )
  
  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_true(is.list(vim$results_all))
})

test_that("varimpact handles single column data frame", {
  # Create test data
  set.seed(1)
  N <- 50
  X_single <- data.frame(x1 = rnorm(N))
  Y <- rbinom(N, 1, plogis(0.2 * X_single$x1))
  
  # Test that single column data frame is processed with warning
  expect_warning(
    vim <- varimpact(Y = Y, data = X_single, verbose = FALSE, V = 2L),
    "Using single variable for variable importance analysis"
  )
  
  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_true(is.list(vim$results_all))
})

test_that("varimpact handles single factor column", {
  # Create test data
  set.seed(1)
  N <- 50
  X_factor <- data.frame(x1 = factor(sample(c("A", "B", "C"), N, replace = TRUE)))
  Y <- rbinom(N, 1, plogis(0.2 * as.numeric(X_factor$x1)))
  
  # Test that single factor column is processed with warning
  expect_warning(
    vim <- varimpact(Y = Y, data = X_factor, verbose = FALSE, V = 2L),
    "Using single variable for variable importance analysis"
  )
  
  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_true(is.list(vim$results_all))
})

test_that("varimpact still works with multiple columns (regression test)", {
  # Create test data
  set.seed(1)
  N <- 50
  X_multi <- data.frame(
    x1 = rnorm(N),
    x2 = rnorm(N),
    x3 = factor(sample(c("A", "B"), N, replace = TRUE))
  )
  Y <- rbinom(N, 1, plogis(0.2 * X_multi$x1 + 0.1 * X_multi$x2))
  
  # Test that multiple columns work without warning
  expect_silent(
    vim <- varimpact(Y = Y, data = X_multi, verbose = FALSE, V = 2L)
  )
  
  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_true(is.list(vim$results_all))
})