library(varimpact)
library(testthat)

context("Recovery when a bin's TMLE estimation fails")

# R/vim-numerics.R and R/vim-factors.R wrap each bin's estimate_tmle2() call in
# try(), then record the failure and carry on. That recovery path was broken:
# `training_estimates[[bin_j]] = NULL` deletes a list element rather than
# storing a NULL placeholder, so the list never grew and the validation code
# indexed past its end. Any single failing bin took down the whole run with
# "subscript out of bounds" - the error reported in issue #20.
#
# These tests force bins to fail so the recovery path is actually exercised.
# Nothing else reaches it reliably: the failure has to come from inside
# estimate_tmle2(), which is hardened enough that ordinary degenerate data no
# longer trips it.

# with_failing() (helper-failure-injection.R) swaps estimate_tmle2() for a
# version that throws on a chosen subset of calls, and restores it afterwards.
with_failing_tmle = function(should_fail, expr) {
  with_failing("estimate_tmle2", should_fail, expr)
}

run_varimpact = function(data, Y, verbose = FALSE, ...) {
  varimpact(Y = Y, data = data, V = 2L,
            Q.library = c("SL.mean", "SL.glm"),
            g.library = c("SL.mean", "SL.glm"),
            verbose = verbose, ...)
}

set.seed(1, "L'Ecuyer-CMRG")
# Large enough, with enough signal, that the happy path yields real results -
# otherwise "it did not crash" would be indistinguishable from "it estimated
# nothing", and the tests below could pass on a broken build.
n = 250
X_num = data.frame(x1 = rnorm(n), x2 = rnorm(n))
Y_bin = rbinom(n, 1, plogis(1.5 * X_num$x1))

X_fac = data.frame(
  f1 = as.factor(sample(c("a", "b", "c"), n, replace = TRUE)),
  f2 = as.factor(sample(c("a", "b", "c"), n, replace = TRUE))
)

future::plan("sequential")

# The case above needs a patched estimate_tmle2(). This one does not: a factor
# level rare enough that some validation folds contain none of it leaves that
# fold with zero rows for that bin, which is the same recovery path reached
# from ordinary data. On master this crashes with "arguments imply differing
# number of rows: 0, 1". It needs V > 2, which is why the package default of
# V = 2 hid it: at V = 2 every fold is large enough to contain the level.
test_that("a rare factor level does not crash multi-fold runs", {
  set.seed(11, "L'Ecuyer-CMRG")
  n = 300
  Y = rbinom(n, 1, 0.5)
  rare = 12
  rest = n - rare
  levels_vec = c(rep("rare", rare),
                 rep("a", rest %/% 2),
                 rep("b", rest - rest %/% 2))
  X = data.frame(f = as.factor(sample(levels_vec)),
                 g = as.factor(sample(c("x", "y"), n, replace = TRUE)))
  vim = varimpact(Y = Y, data = X, V = 5L,
                  Q.library = c("SL.mean", "SL.glm"),
                  g.library = c("SL.mean", "SL.glm"),
                  verbose = FALSE)
  expect_s3_class(vim, "varimpact")
  # Not merely "did not crash": the other levels still have to be estimated.
  expect_true(nrow(vim$results_all) > 0)
})

test_that("numerics survive a failure in some bins", {
  # Every other bin fails, so both the failing and succeeding branches run and
  # the surviving bins still have to line up with their labels.
  vim = with_failing_tmle(function(i) i %% 2 == 0,
                          run_varimpact(X_num, Y_bin))
  expect_s3_class(vim, "varimpact")
})

test_that("numerics survive a failure in the very first bin", {
  # The original bug was worst here: assigning NULL to [[1]] of an empty list
  # left it empty, so the next line indexed element 1 of a zero-length list.
  vim = with_failing_tmle(function(i) i == 1L,
                          run_varimpact(X_num, Y_bin))
  expect_s3_class(vim, "varimpact")
})

test_that("numerics report a failed fold when every bin fails", {
  # All thetas are NA, so which.max() returns integer(0). Comparing or indexing
  # that errored instead of recording the failed fold.
  vim = with_failing_tmle(function(i) TRUE,
                          run_varimpact(X_num, Y_bin))
  expect_s3_class(vim, "varimpact")
  # No variable can be estimated, so nothing survives into the results.
  expect_true(is.null(vim$results_all) || nrow(vim$results_all) == 0)
})

test_that("factors survive a failure in some bins", {
  vim = with_failing_tmle(function(i) i %% 2 == 0,
                          run_varimpact(X_fac, Y_bin))
  expect_s3_class(vim, "varimpact")
})

test_that("factors report a failed fold when every bin fails", {
  vim = with_failing_tmle(function(i) TRUE,
                          run_varimpact(X_fac, Y_bin))
  expect_s3_class(vim, "varimpact")
  expect_true(is.null(vim$results_all) || nrow(vim$results_all) == 0)
})

test_that("the happy path is unaffected and estimate_tmle2 is restored", {
  # Guards against a failed restore leaking into later test files, and confirms
  # the guards added for the failure path did not change ordinary behaviour.
  vim = run_varimpact(X_num, Y_bin)
  expect_s3_class(vim, "varimpact")
  expect_true(nrow(vim$results_all) > 0)
})

test_that("training-bin failures are marked in verbose output", {
  for (X in list(X_num, X_fac)) {
    out = capture.output(
      vim <- with_failing_tmle(function(i) i %% 2 == 0,
                               run_varimpact(X, Y_bin, verbose = TRUE)))
    expect_s3_class(vim, "varimpact")
    # Each failed bin prints an X while the training estimates are running.
    expect_true(any(grepl("X", out, fixed = TRUE)))
  }
  # When every bin fails, compile_results() explains what it found.
  out = capture.output(suppressWarnings(
    vim <- with_failing_tmle(function(i) TRUE,
                             run_varimpact(X_num, Y_bin, verbose = TRUE))))
  expect_true(any(grepl("EY1 lengths:", out, fixed = TRUE)))
})
