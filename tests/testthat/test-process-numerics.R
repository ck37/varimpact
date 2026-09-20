library(varimpact)
library(testthat)

# process_numerics() is only ever called by varimpact() with continuous
# variables that have more distinct values than bins. Called directly, it
# also has to handle a variable with few distinct values (binned on its own
# values rather than quantiles) and report the columns it drops.

context("process_numerics")

test_that("a variable with no more distinct values than bins keeps them as its bins", {
  set.seed(31, "L'Ecuyer-CMRG")
  n = 100
  data = data.frame(few = sample(0:2, n, replace = TRUE), many = rnorm(n))
  out = varimpact:::process_numerics(data, quantile_probs_numeric = c(0.1, 0.9),
                                     miss.cut = 0.5, bins_numeric = 10L,
                                     impute = "median")
  expect_equal(out$num_numeric, 2L)
  # Three distinct values, three bins; the continuous variable gets ten.
  # data.cont.dist holds the bin index of every observation.
  bins = function(col) length(unique(na.omit(out$data.cont.dist[, col])))
  expect_equal(bins("few"), 3L)
  expect_equal(bins("many"), 10L)
})

test_that("verbose reports the numerics dropped for lack of variation", {
  set.seed(32, "L'Ecuyer-CMRG")
  n = 100
  data = data.frame(flat = rep(1, n), many = rnorm(n))
  expect_output(
    out <- varimpact:::process_numerics(data, quantile_probs_numeric = c(0.1, 0.9),
                                        miss.cut = 0.5, bins_numeric = 5L,
                                        impute = "median", verbose = TRUE),
    "Dropped 1 numerics due to lack of variation")
  expect_equal(out$num_numeric, 1L)
})
