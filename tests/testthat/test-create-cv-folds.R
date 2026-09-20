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

test_that("assignment interleaves by row order, as the cvTools version did", {
  # create_cv_folds() read cvTools::cvFolds(type = "random")$which, which is
  # rep(seq_len(V), length.out = n): the fold of the permuted observation, not
  # of the original one. The permutation was never used, so the assignment was
  # deterministic all along. This pins that so the cvTools removal changes no
  # result; see assign_folds() for why it is not randomized here.
  Y = c(rep(0, 7), rep(1, 5))
  folds = varimpact:::create_cv_folds(V = 3L, Y = Y)
  expect_identical(folds[Y == 0], rep_len(1:3, 7))
  expect_identical(folds[Y == 1], rep_len(1:3, 5))

  set.seed(1)
  expect_identical(varimpact:::create_cv_folds(V = 2L, Y = Y),
                   varimpact:::create_cv_folds(V = 2L, Y = Y))
  expect_identical(varimpact:::assign_folds(10L, 4L), rep_len(1:4, 10))
})

test_that("the RNG advances as it did with cvTools, so seeded fits are unchanged", {
  # cvFolds() drew an unused permutation of each stratum. assign_folds() makes
  # the same draw so that whatever draws from the RNG afterward gets the same
  # numbers as before, and a fit from a given seed reproduces the previous
  # version's. Every stratum is one sample.int(n) call, and the RNG state
  # after create_cv_folds() must equal the state after those calls.
  Y = c(rep(0, 7), rep(1, 5), NA, NA)
  set.seed(3, "L'Ecuyer-CMRG")
  invisible(varimpact:::create_cv_folds(V = 2L, Y = Y))
  after_folds = runif(1)

  set.seed(3, "L'Ecuyer-CMRG")
  for (n in c(7L, 5L, 2L)) sample.int(n)
  expect_identical(runif(1), after_folds)
})

test_that("a stratum smaller than V does not error", {
  # cvTools::cvFolds() refused K > n. A rare outcome level with fewer
  # observations than folds is still spread over the folds it can reach.
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
