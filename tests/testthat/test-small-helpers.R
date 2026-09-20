library(varimpact)
library(testthat)

# The remaining branches of the small helpers: the no-results print, the
# quantile helpers' warning and verbose paths.

context("small helpers")

test_that("print.varimpact says so when there are no results at all", {
  empty = structure(list(results_consistent = NULL, results_all = NULL),
                    class = "varimpact")
  expect_output(print(empty), "No results could be calculated")
})

test_that("print.varimpact falls back to all results when none are consistent", {
  vim = structure(list(results_consistent = data.frame(),
                       results_all = data.frame(Estimate = 0.1, row.names = "x1")),
                  class = "varimpact")
  expect_output(print(vim), "No significant and consistent results")
  expect_output(print(vim), "x1")
})

test_that("quantiles_equivalent warns when given the wrong number of probabilities", {
  expect_warning(quantiles_equivalent(1:10, quantile_probs = c(0.1, 0.5, 0.9)),
                 "2-element vector")
})

test_that("quantiles_equivalent accepts a factor", {
  expect_true(quantiles_equivalent(factor(rep("a", 20))))
  expect_false(quantiles_equivalent(factor(rep(c("a", "b"), 10))))
})

test_that("restrict_by_quantiles names the columns it drops when verbose", {
  data = data.frame(keep = rnorm(50), flat = rep(1, 50))
  expect_output(out <- restrict_by_quantiles(data, verbose = TRUE), "dropping flat")
  expect_identical(names(out), "keep")
  # And is silent when there is nothing to drop.
  expect_silent(restrict_by_quantiles(data["keep"], verbose = TRUE))
})
