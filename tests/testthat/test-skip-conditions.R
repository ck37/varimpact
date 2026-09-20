library(varimpact)
library(testthat)

# vim_factors() and vim_numerics() skip a variable in a fold when it cannot be
# estimated: too few events for a binary outcome, a factor with one level left
# in the training fold, or HOPACH reducing the adjustment set to constants.
# Each records a message and moves on. These reach each of those branches
# and the "nothing could be estimated at all" outcome in compile_results().

context("Skip conditions in the per-fold estimation")

future::plan("sequential")
lib = c("SL.mean", "SL.glm")

run_verbose = function(...) {
  out = capture.output(vim <- varimpact(V = 2L, verbose = TRUE, Q.library = lib,
                                        g.library = lib, ...))
  list(vim = vim, text = paste(out, collapse = "\n"))
}

test_that("a binary outcome with too few events skips the variable (minYs)", {
  set.seed(21, "L'Ecuyer-CMRG")
  n = 100
  X = data.frame(f = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  # About 12 events in total, so roughly 6 per training fold, under the
  # default minYs = 15.
  Y = rbinom(n, 1, 0.12)
  res = run_verbose(Y = Y, data = X)
  expect_s3_class(res$vim, "varimpact")
  expect_true(grepl("due to minY constraint", res$text, fixed = TRUE))
})

test_that("a factor with a single level left in the training fold is skipped", {
  set.seed(22, "L'Ecuyer-CMRG")
  n = 100
  # One observation of "b". Whichever fold it lands in, the other fold's
  # training data has only "a", so that fold cannot estimate anything.
  # quantile_probs_factor = NULL keeps process_factors() from dropping the
  # variable for lack of variation before it gets that far.
  X = data.frame(f = factor(c("b", rep("a", n - 1))),
                 g = factor(sample(c("x", "y"), n, replace = TRUE)))
  Y = rbinom(n, 1, 0.5)
  res = run_verbose(Y = Y, data = X, quantile_probs_factor = NULL)
  expect_s3_class(res$vim, "varimpact")
  expect_true(grepl("due to lack of variation", res$text, fixed = TRUE))
})

test_that("an adjustment set reduced to constants is skipped, for numerics and factors", {
  # reduce_dimensions() never returns constant columns itself, so make it.
  make_constant = function(original) {
    function(data, newX = NULL, max_variables, verbose = FALSE) {
      out = original(data, newX, max_variables, verbose)
      out$data[] = 1
      out$newX[] = 1
      out
    }
  }
  set.seed(23, "L'Ecuyer-CMRG")
  n = 120
  X_num = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  X_fac = data.frame(f1 = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
                     f2 = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  Y = rbinom(n, 1, 0.5)

  res = with_replaced("reduce_dimensions", make_constant, run_verbose(Y = Y, data = X_num))
  expect_s3_class(res$vim, "varimpact")
  expect_true(grepl("HOPACH reduced W to all constant columns", res$text, fixed = TRUE))
  expect_true(grepl("Constant columns (", res$text, fixed = TRUE))

  res = with_replaced("reduce_dimensions", make_constant, run_verbose(Y = Y, data = X_fac))
  expect_s3_class(res$vim, "varimpact")
  expect_true(grepl("HOPACH reduced W to all constant columns", res$text, fixed = TRUE))
})

test_that("data with no usable variables yields an empty result and a warning", {
  set.seed(24, "L'Ecuyer-CMRG")
  n = 80
  Y = rbinom(n, 1, 0.5)
  constant = data.frame(a = rep(1, n), b = rep(2, n))
  out = capture.output(expect_warning(
    vim <- varimpact(Y = Y, data = constant, V = 2L, verbose = TRUE,
                     Q.library = "SL.mean", g.library = "SL.mean"),
    "No variable importance estimates could be calculated"))
  expect_s3_class(vim, "varimpact")
  expect_null(vim$results_all)
  expect_true(any(grepl("No variable importance estimates", out, fixed = TRUE)))
})
