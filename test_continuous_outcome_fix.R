#!/usr/bin/env Rscript

# Test script to demonstrate the fix for continuous outcome estimates
# Issue #8: Transform continuous outcome estimates back to original scale

library(varimpact)

# Create test dataset with continuous outcome
set.seed(1, "L'Ecuyer-CMRG")
N <- 100
num_normal <- 3
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))

# Create a continuous outcome with a known range
# Let's make Y have values roughly between 10 and 50
Y_continuous <- 30 + 10 * (.2*X[, 1] + .1*X[, 2] - .2*X[, 3]) + rnorm(N, 0, 2)

cat("Original Y range:", range(Y_continuous), "\n")
cat("Original Y mean:", mean(Y_continuous), "\n")

# Run varimpact with gaussian family
vim <- varimpact(Y = Y_continuous, data = X, 
                 family = "gaussian",
                 V = 2L,
                 verbose = TRUE,
                 Q.library = c("SL.mean", "SL.glm"),
                 g.library = c("SL.mean", "SL.glm"))

cat("\nResults:\n")
print(vim$results_all)

# Check if the estimates are on the original scale
estimates <- vim$results_all$Estimate
cat("\nEstimate range:", range(estimates, na.rm = TRUE), "\n")
cat("Estimate mean:", mean(estimates, na.rm = TRUE), "\n")

# The estimates should be on the original scale (around 10-50 range)
# not on the [0, 1] scale
if (all(estimates >= 0 & estimates <= 1, na.rm = TRUE)) {
  cat("ERROR: Estimates appear to be on [0, 1] scale, not original scale!\n")
} else {
  cat("SUCCESS: Estimates appear to be on the original scale!\n")
}

# Also check the individual fold estimates
cat("\nFold estimates:\n")
if (!is.null(vim$results_by_fold)) {
  fold_estimates <- as.matrix(vim$results_by_fold[, grep("Est_v", colnames(vim$results_by_fold))])
  cat("Fold estimate range:", range(fold_estimates, na.rm = TRUE), "\n")
  
  if (all(fold_estimates >= 0 & fold_estimates <= 1, na.rm = TRUE)) {
    cat("ERROR: Fold estimates appear to be on [0, 1] scale!\n")
  } else {
    cat("SUCCESS: Fold estimates appear to be on the original scale!\n")
  }
}