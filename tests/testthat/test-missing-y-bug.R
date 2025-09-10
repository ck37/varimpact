# Test for GitHub issue #32: Bug with missing Y values
# https://github.com/ck37/varimpact/issues/32
#
# The bug is in vim-factors.R lines 181-186:
# if (sum(deltat == 0) < 10) {
#   Yt = Yt[deltat == 1]
#   At = At[deltat == 1]
#   Wtsht = Wtsht[deltat == 1, , drop = FALSE]
#   deltat = deltat[deltat == 1]
# }
#
# This code only cleans up missing values when there are <10 missing observations,
# but when there are >=10 missing values, the cleanup doesn't happen and downstream
# TMLE code fails. The entire block should be removed per the TODO comment.

library(testthat)
library(varimpact)

context("Missing Y values bug (Issue #32)")

test_that("varimpact fails with >10 missing Y values due to inconsistent missing data handling", {
  # Reproduce the exact simulation from the GitHub issue
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 200
  num_normal <- 4
  X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
  Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
  
  # Add some missing data to X so we can test imputation
  for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA
  
  # Add missing data to Y - exactly as in the GitHub issue
  # This creates 11 missing Y values, which triggers the bug
  Y[c(4,6,7,8,11,15,20,21,28,32,72)] <- NA
  
  # Verify we have exactly 11 missing Y values (>10)
  missing_y_count <- sum(is.na(Y))
  expect_equal(missing_y_count, 11, info = paste("Expected 11 missing Y values, got", missing_y_count))
  
  # This should fail with the current buggy code because:
  # 1. deltat will have 11 zeros (missing observations)
  # 2. sum(deltat == 0) = 11, which is NOT < 10
  # 3. So the cleanup code doesn't run
  # 4. But TMLE downstream expects clean data and crashes
  expect_error({
    vim <- varimpact(Y = Y, data = X, 
                     V = 2L,  # Use fewer folds for faster testing
                     Q.library = c("SL.mean", "SL.glm"),
                     g.library = c("SL.mean", "SL.glm"),
                     verbose = FALSE)
  }, info = "Current code should fail with >10 missing Y values")
})

test_that("varimpact works with <10 missing Y values (but uses problematic code path)", {
  # Test the case where the problematic code path is taken (< 10 missing)
  # This demonstrates the inconsistent behavior: works with <10 missing but fails with >=10
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 200
  num_normal <- 4
  X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
  Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
  
  # Add some missing data to X
  for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA
  
  # Add only a few missing Y values (< 10)
  Y[c(4,6,7)] <- NA
  
  # Verify we have exactly 3 missing Y values (<10)
  missing_y_count <- sum(is.na(Y))
  expect_equal(missing_y_count, 3, info = paste("Expected 3 missing Y values, got", missing_y_count))
  
  # This should work because sum(deltat == 0) = 3 < 10, so cleanup code runs
  # But this demonstrates the inconsistent behavior that should be fixed
  expect_silent({
    vim <- varimpact(Y = Y, data = X, 
                     V = 2L,  # Use fewer folds for faster testing
                     Q.library = c("SL.mean", "SL.glm"),
                     g.library = c("SL.mean", "SL.glm"),
                     verbose = FALSE)
  })
})

test_that("edge case: exactly 10 missing Y values", {
  # Test the exact boundary condition
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 200
  num_normal <- 4
  X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
  Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
  
  # Add some missing data to X
  for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA
  
  # Add exactly 10 missing Y values - this should trigger the bug
  # because sum(deltat == 0) = 10, which is NOT < 10
  Y[c(4,6,7,8,11,15,20,21,28,32)] <- NA
  
  # Verify we have exactly 10 missing Y values
  missing_y_count <- sum(is.na(Y))
  expect_equal(missing_y_count, 10, info = paste("Expected 10 missing Y values, got", missing_y_count))
  
  # This should fail because sum(deltat == 0) = 10 is NOT < 10
  # So cleanup code doesn't run, but TMLE expects clean data
  expect_error({
    vim <- varimpact(Y = Y, data = X, 
                     V = 2L,  # Use fewer folds for faster testing
                     Q.library = c("SL.mean", "SL.glm"),
                     g.library = c("SL.mean", "SL.glm"),
                     verbose = FALSE)
  }, info = "Should fail with exactly 10 missing Y values due to >= 10 condition")
})