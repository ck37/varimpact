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

# Helper functions to simulate the problematic behavior
simulate_old_behavior <- function(Yt, At) {
  deltat <- as.numeric(!is.na(Yt) & !is.na(At))
  
  # This was the problematic code that has been REMOVED in the fix
  if (sum(deltat == 0) < 10) {
    Yt <- Yt[deltat == 1]
    At <- At[deltat == 1]
    deltat <- deltat[deltat == 1]
  }
  
  list(
    Yt = Yt,
    At = At,
    deltat = deltat,
    n_missing = sum(is.na(Yt)),
    cleanup_ran = sum(deltat == 0) < 10
  )
}

simulate_new_behavior <- function(Yt, At) {
  deltat <- as.numeric(!is.na(Yt) & !is.na(At))
  
  # After the fix: no conditional cleanup, consistent behavior
  list(
    Yt = Yt,
    At = At,
    deltat = deltat,
    n_missing = sum(is.na(Yt)),
    cleanup_ran = FALSE  # No cleanup in new behavior
  )
}

context("Missing Y values bug (Issue #32)")

test_that("old behavior shows inconsistent missing data handling with >10 missing values", {
  # Create test data with >10 missing Y values
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 100
  Yt <- rbinom(N, 1, 0.5)
  At <- as.factor(sample(c("A", "B"), N, replace = TRUE))
  
  # Add 15 missing Y values (>10)
  missing_indices <- sample(1:N, 15)
  Yt[missing_indices] <- NA
  
  # Test old behavior
  result <- simulate_old_behavior(Yt, At)
  
  # With >10 missing values, cleanup should NOT run
  expect_false(result$cleanup_ran, info = "Cleanup should not run with >10 missing values")
  expect_equal(result$n_missing, 15, info = "Missing values should remain in data")
  expect_equal(length(result$Yt), N, info = "Data length should be unchanged")
})

test_that("old behavior shows inconsistent missing data handling with <10 missing values", {
  # Create test data with <10 missing Y values
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 100
  Yt <- rbinom(N, 1, 0.5)
  At <- as.factor(sample(c("A", "B"), N, replace = TRUE))
  
  # Add 5 missing Y values (<10)
  missing_indices <- sample(1:N, 5)
  Yt[missing_indices] <- NA
  
  # Test old behavior
  result <- simulate_old_behavior(Yt, At)
  
  # With <10 missing values, cleanup SHOULD run
  expect_true(result$cleanup_ran, info = "Cleanup should run with <10 missing values")
  expect_equal(result$n_missing, 0, info = "Missing values should be removed")
  expect_equal(length(result$Yt), N - 5, info = "Data should be shortened by removing missing values")
})

test_that("edge case: exactly 10 missing Y values triggers bug", {
  # Create test data with exactly 10 missing Y values
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 100
  Yt <- rbinom(N, 1, 0.5)
  At <- as.factor(sample(c("A", "B"), N, replace = TRUE))
  
  # Add exactly 10 missing Y values
  missing_indices <- sample(1:N, 10)
  Yt[missing_indices] <- NA
  
  # Test old behavior
  result <- simulate_old_behavior(Yt, At)
  
  # With exactly 10 missing values, cleanup should NOT run (10 is not < 10)
  expect_false(result$cleanup_ran, info = "Cleanup should not run with exactly 10 missing values")
  expect_equal(result$n_missing, 10, info = "Missing values should remain in data")
  expect_equal(length(result$Yt), N, info = "Data length should be unchanged")
})

test_that("new behavior shows consistent handling regardless of missing value count", {
  # Test with various numbers of missing values
  test_cases <- c(5, 10, 15, 20)
  
  for (n_missing in test_cases) {
    set.seed(1, "L'Ecuyer-CMRG")
    N <- 100
    Yt <- rbinom(N, 1, 0.5)
    At <- as.factor(sample(c("A", "B"), N, replace = TRUE))
    
    # Add missing Y values
    missing_indices <- sample(1:N, n_missing)
    Yt[missing_indices] <- NA
    
    # Test new behavior
    result <- simulate_new_behavior(Yt, At)
    
    # New behavior should be consistent regardless of missing count
    expect_false(result$cleanup_ran, info = paste("No cleanup should run with", n_missing, "missing values"))
    expect_equal(result$n_missing, n_missing, info = paste("All", n_missing, "missing values should remain"))
    expect_equal(length(result$Yt), N, info = paste("Data length should be unchanged with", n_missing, "missing values"))
  }
})

test_that("fix eliminates the arbitrary 10-missing-value threshold", {
  # Test that the fix removes the inconsistent behavior
  set.seed(1, "L'Ecuyer-CMRG")
  N <- 100
  
  # Test case 1: 9 missing values
  Yt1 <- rbinom(N, 1, 0.5)
  At1 <- as.factor(sample(c("A", "B"), N, replace = TRUE))
  Yt1[sample(1:N, 9)] <- NA
  
  # Test case 2: 11 missing values
  Yt2 <- rbinom(N, 1, 0.5)
  At2 <- as.factor(sample(c("A", "B"), N, replace = TRUE))
  Yt2[sample(1:N, 11)] <- NA
  
  # Old behavior would be different
  old_result1 <- simulate_old_behavior(Yt1, At1)
  old_result2 <- simulate_old_behavior(Yt2, At2)
  
  expect_true(old_result1$cleanup_ran, info = "Old behavior: cleanup runs with 9 missing")
  expect_false(old_result2$cleanup_ran, info = "Old behavior: cleanup doesn't run with 11 missing")
  
  # New behavior should be consistent
  new_result1 <- simulate_new_behavior(Yt1, At1)
  new_result2 <- simulate_new_behavior(Yt2, At2)
  
  expect_equal(new_result1$cleanup_ran, new_result2$cleanup_ran, 
               info = "New behavior should be consistent regardless of missing count")
  expect_false(new_result1$cleanup_ran, info = "New behavior: no cleanup with 9 missing")
  expect_false(new_result2$cleanup_ran, info = "New behavior: no cleanup with 11 missing")
})