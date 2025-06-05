library(varimpact)
library(testthat)

context("Continuous outcome scale transformation")

test_that("continuous outcomes are reported on original scale, not [0,1] scale", {
  # Set multicore-compatible seed
  set.seed(42, "L'Ecuyer-CMRG")
  
  # Create test dataset with continuous outcome in a known range
  N <- 100
  num_vars <- 3
  X <- as.data.frame(matrix(rnorm(N * num_vars), N, num_vars))
  colnames(X) <- paste0("X", 1:num_vars)
  
  # Create a continuous outcome with values roughly between 20 and 80
  # This gives us a clear range that's far from [0, 1]
  Y_continuous <- 50 + 15 * (.3*X[, 1] + .2*X[, 2] - .1*X[, 3]) + rnorm(N, 0, 3)
  
  # Verify our outcome is in the expected range
  expect_true(min(Y_continuous) > 10, "Y should be well above 0")
  expect_true(max(Y_continuous) < 100, "Y should be well below 100")
  expect_true(max(Y_continuous) - min(Y_continuous) > 10, "Y should have substantial range")
  
  # Run varimpact with gaussian family
  # Use sequential execution and simple libraries for faster testing
  future::plan("sequential")
  
  vim <- varimpact(Y = Y_continuous, 
                   data = X, 
                   family = "gaussian",
                   V = 2L,  # Use fewer folds for faster testing
                   verbose = FALSE,
                   Q.library = c("SL.mean", "SL.glm"),
                   g.library = c("SL.mean", "SL.glm"),
                   bins_numeric = 2L)  # Fewer bins for faster testing
  
  # Check that the function completed successfully
  expect_s3_class(vim, "varimpact")
  expect_true(!is.null(vim$results_all))
  
  # Extract estimates
  estimates <- vim$results_all$Estimate
  estimates <- estimates[!is.na(estimates)]
  
  # The key test: estimates should NOT be on [0, 1] scale
  # If the bug exists, all estimates would be between 0 and 1
  # With the fix, estimates should be on the original scale
  
  # Check that not all estimates are in [0, 1] range
  # (some estimates might legitimately be in [0, 1] by chance, but not all)
  estimates_in_01_range <- estimates >= 0 & estimates <= 1
  proportion_in_01 <- mean(estimates_in_01_range)
  
  # If more than 90% of estimates are in [0, 1], likely the bug still exists
  expect_lt(proportion_in_01, 0.9, 
            paste("Too many estimates in [0,1] range. Estimates:", 
                  paste(round(estimates, 3), collapse = ", ")))
  
  # Additional check: the range of estimates should be reasonable
  # relative to the original outcome range
  original_range <- max(Y_continuous) - min(Y_continuous)
  estimate_range <- max(estimates) - min(estimates)
  
  # The estimate range should be a reasonable fraction of the original range
  # (not tiny like it would be if stuck on [0, 1] scale)
  expect_gt(estimate_range, original_range * 0.01,
            "Estimate range seems too small relative to original outcome range")
  
  # Test fold-specific estimates if available
  if (!is.null(vim$results_by_fold)) {
    fold_est_cols <- grep("Est_v", colnames(vim$results_by_fold), value = TRUE)
    if (length(fold_est_cols) > 0) {
      fold_estimates <- as.matrix(vim$results_by_fold[, fold_est_cols])
      fold_estimates <- fold_estimates[!is.na(fold_estimates)]
      
      if (length(fold_estimates) > 0) {
        fold_estimates_in_01 <- fold_estimates >= 0 & fold_estimates <= 1
        fold_proportion_in_01 <- mean(fold_estimates_in_01)
        
        expect_lt(fold_proportion_in_01, 0.9,
                  "Too many fold estimates in [0,1] range")
      }
    }
  }
})

test_that("binary outcomes still work correctly", {
  # Set seed for reproducibility
  set.seed(43, "L'Ecuyer-CMRG")
  
  # Create test dataset with binary outcome
  N <- 100
  X <- as.data.frame(matrix(rnorm(N * 2), N, 2))
  colnames(X) <- c("X1", "X2")
  
  # Create binary outcome
  Y_binary <- rbinom(N, 1, plogis(.5*X[, 1] - .3*X[, 2]))
  
  # Run varimpact with binomial family
  future::plan("sequential")
  
  vim <- varimpact(Y = Y_binary, 
                   data = X, 
                   family = "binomial",
                   V = 2L,
                   verbose = FALSE,
                   Q.library = c("SL.mean", "SL.glm"),
                   g.library = c("SL.mean", "SL.glm"))
  
  # Check that the function completed successfully
  expect_s3_class(vim, "varimpact")
  expect_true(!is.null(vim$results_all))
  
  # For binary outcomes, estimates should be in a reasonable range
  # (typically between -1 and 1 for risk differences)
  estimates <- vim$results_all$Estimate
  estimates <- estimates[!is.na(estimates)]
  
  if (length(estimates) > 0) {
    expect_true(all(estimates >= -1 & estimates <= 1),
                "Binary outcome estimates should be reasonable risk differences")
  }
})