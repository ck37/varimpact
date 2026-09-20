library(varimpact)
library(testthat)

context("create_cv_folds")

# Fold sizes within a stratum may differ by at most one.
expect_balanced = function(folds, V) {
  sizes = tabulate(folds, nbins = V)
  expect_lte(max(sizes) - min(sizes), 1L)
}

test_that("a binary outcome is split evenly within each level", {
  set.seed(1, "L'Ecuyer-CMRG")
  Y = c(rep(0, 45), rep(1, 25))
  folds = varimpact:::create_cv_folds(V = 3L, Y = Y)

  expect_equal(length(folds), length(Y))
  expect_false(anyNA(folds))
  expect_setequal(unique(folds), 1:3)
  expect_balanced(folds[Y == 0], 3L)
  expect_balanced(folds[Y == 1], 3L)
})

test_that("a non-binary outcome is split evenly without stratification", {
  set.seed(2, "L'Ecuyer-CMRG")
  Y = rnorm(41)
  folds = varimpact:::create_cv_folds(V = 4L, Y = Y)

  expect_equal(length(folds), length(Y))
  expect_setequal(unique(folds), 1:4)
  expect_balanced(folds, 4L)
})

test_that("assignment is random, and reproducible from the seed", {
  Y = rbinom(60, 1, 0.5)
  set.seed(3, "L'Ecuyer-CMRG")
  first = varimpact:::create_cv_folds(V = 2L, Y = Y)
  set.seed(3, "L'Ecuyer-CMRG")
  again = varimpact:::create_cv_folds(V = 2L, Y = Y)
  set.seed(4, "L'Ecuyer-CMRG")
  other = varimpact:::create_cv_folds(V = 2L, Y = Y)

  expect_identical(first, again)
  expect_false(identical(first, other))
  # Not the interleaved pattern 1, 2, 1, 2, ... that the cvTools-based version
  # produced, which never used the permutation it drew.
  expect_false(identical(first[Y == 0], rep_len(1:2, sum(Y == 0))))
})

test_that("a stratum smaller than V does not error", {
  # cvTools::cvFolds() refused K > n. A rare outcome level with fewer
  # observations than folds should still get spread over the folds it can.
  set.seed(5, "L'Ecuyer-CMRG")
  Y = c(rep(0, 30), 1, 1)
  folds = varimpact:::create_cv_folds(V = 5L, Y = Y)
  expect_false(anyNA(folds))
  expect_equal(length(unique(folds[Y == 1])), 2L)
  expect_balanced(folds[Y == 0], 5L)

  expect_equal(length(varimpact:::assign_folds(1L, 3L)), 1L)
})

test_that("verbose prints the fold breakdown", {
  set.seed(6, "L'Ecuyer-CMRG")
  Y = rbinom(40, 1, 0.5)
  expect_output(varimpact:::create_cv_folds(V = 2L, Y = Y, verbose = TRUE),
                "Cross-validation fold breakdown")
})
