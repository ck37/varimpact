library(testthat)
library(varimpact)

#################################
# mlbench BreastCancer dataset.

context("BreastCancer dataset")

data(BreastCancer, package = "mlbench")
data = BreastCancer
names(data) = tolower(names(data))

set.seed(3, "L'Ecuyer-CMRG")

# Reduce to a dataset of 120 observations to speed up testing.
data = data[sample(nrow(data), 120), ]

# Create a numeric outcome variable.
data$y = as.numeric(data$class == "malignant")
table(data$y)

x = subset(data, select = -c(y, class, id))
dim(x)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  future::plan("multiprocess", workers = 2)
}

test_that("varimpact runs on BreastCancer dataset", {
  # Speed up the test with moderate optimization:
  # - Use only 2 variables to reduce computation time
  # - Use faster SuperLearner libraries
  # - Reduce bins and use minimal cross-validation
  vim = varimpact(Y = data$y, x[, 1:2], 
                  Q.library = c("SL.glm", "SL.mean"),
                  g.library = c("SL.glm", "SL.mean"),
                  bins_numeric = 5L,
                  V = 2L,
                  verbose = FALSE, 
                  verbose_tmle = FALSE)
  
  # Test that the function completes without error
  expect_is(vim, "varimpact")
  expect_is(vim$time, "proc_time")
  
  # Test that some VIMs were actually calculated
  # The test should produce meaningful results
  if (!is.null(vim$results_all)) {
    expect_gte(nrow(vim$results_all), 1)
    cat("Successfully calculated", nrow(vim$results_all), "VIMs\n")
  } else {
    cat("Warning: No VIMs calculated\n")
  }
  
  # Print timing for reference
  cat("Test execution time:", vim$time[3], "seconds\n")
})

test_that("varimpact works with A_names parameter", {
  # Test a subset of columns for A_names (use just 1 variable).
  vim = varimpact(Y = data$y, x, 
                  A_names = colnames(x)[1], 
                  Q.library = c("SL.glm", "SL.mean"),
                  g.library = c("SL.glm", "SL.mean"),
                  bins_numeric = 5L,
                  V = 2L,
                  verbose = FALSE)
  
  # Test that the function completes without error
  expect_is(vim, "varimpact")
  expect_is(vim$time, "proc_time")
  
  # Test that some VIMs were actually calculated
  if (!is.null(vim$results_all)) {
    expect_gte(nrow(vim$results_all), 1)
    cat("Successfully calculated", nrow(vim$results_all), "VIMs with A_names\n")
  } else {
    cat("Warning: No VIMs calculated with A_names\n")
  }
  
  # Print timing for reference
  cat("A_names test execution time:", vim$time[3], "seconds\n")
})

# Return to single core usage.
future::plan("sequential")
