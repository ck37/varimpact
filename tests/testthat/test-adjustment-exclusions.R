# Tests for the adjustment_exclusions argument of varimpact(), contributed in
# https://github.com/ck37/varimpact/pull/25.
#
# These assert the behaviour that matters - the excluded variables are absent
# from the adjustment set and everything else is still present - rather than
# comparing W_names against a hard-coded vector, which would break whenever the
# missingness-indicator or dimension-reduction machinery changes.

library(testthat)
library(varimpact)

context("Adjustment set exclusions")

# The columns of an adjustment matrix that belong to one input variable.
# Mirrors how varimpact expands data: "V1", "F1XXb"/"F1XXc", "Imiss_V1".
test_that("adjustment_columns() resolves numerics, factors and missingness", {
  w <- c("V1", "V2", "F1XXb", "F1XXc", "Imiss_V1", "Imiss_F1", "F10XXa")

  # A numeric takes its own column and its missingness indicator.
  expect_equal(varimpact:::adjustment_columns("V1", w), c("V1", "Imiss_V1"))

  # A factor takes every indicator column, not a column of its own name.
  expect_equal(varimpact:::adjustment_columns("F1", w),
               c("F1XXb", "F1XXc", "Imiss_F1"))

  # A variable that contributes nothing matches nothing.
  expect_equal(varimpact:::adjustment_columns("V9", w), character(0))
})

test_that("exclude_adjustment_vars() drops only the requested variables", {
  W <- data.frame(V1 = 1:3, V2 = 1:3, F1XXb = 1:3, F1XXc = 1:3, Imiss_V1 = 1:3)

  # Nothing requested for this variable: W is untouched.
  expect_equal(colnames(varimpact:::exclude_adjustment_vars(W, "V1", list())),
               colnames(W))
  expect_equal(
    colnames(varimpact:::exclude_adjustment_vars(W, "V1", list(V2 = "V1"))),
    colnames(W))

  # Excluding a numeric also removes its missingness indicator.
  expect_equal(
    colnames(varimpact:::exclude_adjustment_vars(W, "V2", list(V2 = "V1"))),
    c("V2", "F1XXb", "F1XXc"))

  # Excluding a factor removes every indicator column.
  expect_equal(
    colnames(varimpact:::exclude_adjustment_vars(W, "V1", list(V1 = "F1"))),
    c("V1", "V2", "Imiss_V1"))

  # A data frame is returned even when a single column survives.
  one <- varimpact:::exclude_adjustment_vars(W, "V1", list(V1 = c("V2", "F1")))
  expect_s3_class(one, "data.frame")
  expect_equal(colnames(one), c("V1", "Imiss_V1"))
})

test_that("check_adjustment_exclusions() warns about names not in the data", {
  data_names <- c("V1", "V2", "F1")

  expect_silent(
    varimpact:::check_adjustment_exclusions(list(V1 = "V2"), data_names))
  expect_silent(varimpact:::check_adjustment_exclusions(list(), data_names))

  # The target variable is unknown.
  expect_warning(
    varimpact:::check_adjustment_exclusions(list(typo = "V2"), data_names),
    "not columns of data")

  # The excluded variable is unknown - the dangerous case, because excluding
  # nothing silently would look identical to a successful exclusion.
  expect_warning(
    varimpact:::check_adjustment_exclusions(list(V1 = "V9"), data_names),
    "not be excluded")

  expect_error(
    varimpact:::check_adjustment_exclusions(list("V2"), data_names),
    "named list")
})

test_that("varimpact() honours adjustment_exclusions for numerics and factors", {
  set.seed(1, "L'Ecuyer-CMRG")
  future::plan("sequential")

  N <- 200
  X <- data.frame(V1 = rnorm(N), V2 = rnorm(N), V3 = rnorm(N),
                  F1 = factor(sample(c("a", "b", "c"), N, replace = TRUE)))
  Y <- rbinom(N, 1, plogis(0.3 * X$V1 + 0.2 * X$V2))

  vim <- suppressWarnings(
    varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"),
              adjustment_exclusions = list(V1 = "V2", F1 = c("V1", "V2"))))

  adjustment_set <- function(name) {
    vim$all_vims[[name]]$fold_results[[1]]$bin_results[[1]]$W_names
  }

  # A numeric candidate: V2 is gone, the others remain.
  v1 <- adjustment_set("V1")
  expect_false("V2" %in% v1)
  expect_true("V3" %in% v1)

  # A factor candidate. This is what PR #25 could not do: adjustment_exclusions
  # never reached vim_factors, so the request was silently ignored.
  f1 <- adjustment_set("F1")
  expect_false(any(c("V1", "V2") %in% f1))
  expect_true("V3" %in% f1)

  # A variable with no exclusions keeps its full adjustment set. Use V2 rather
  # than V3: penalized histogramming can collapse a variable to a single level
  # in a fold, in which case it has no bins and records no adjustment set, and
  # that has nothing to do with exclusions.
  v2 <- adjustment_set("V2")
  expect_true(all(c("V1", "V3") %in% v2))
})
