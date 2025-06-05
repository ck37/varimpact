library(testthat)
library(varimpact)

#################################
# mlbench BreastCancer dataset.

context("BreastCancer dataset")

data(BreastCancer, package = "mlbench")
data = BreastCancer
names(data) = tolower(names(data))

set.seed(3, "L'Ecuyer-CMRG")

# Reduce to a dataset of 50 observations to speed up testing.
data = data[sample(nrow(data), 50), ]

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
  # Speed up the test with very aggressive optimization:
  # - Use only 1 variable to minimize computation
  # - Use only SL.mean (fastest possible)
  # - Minimal bins and very low thresholds
  # - Single cross-validation fold
  vim = varimpact(Y = data$y, x[, 1:2], 
                  Q.library = c("SL.mean"), 
                  g.library = c("SL.mean"),
                  bins_numeric = 3L,
                  minCell = 0L,
                  minYs = 1L,
                  V = 2L,
                  verbose = FALSE, 
                  verbose_tmle = FALSE)
  
  # Test that the function completes without error
  expect_is(vim, "varimpact")
  expect_is(vim$time, "proc_time")
  
  # Print timing for reference
  cat("Test execution time:", vim$time[3], "seconds\n")
})

# Return to single core usage.
future::plan("sequential")
