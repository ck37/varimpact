library(varimpact)
library(testthat)

# vim_numerics() and vim_factors() call apply_tmle_to_validation() twice per
# fold - once for the low level, once for the high level - inside try(), and
# record a message when it fails. Those two branches are separate code and
# neither is reachable from ordinary data, so the failures are injected.

context("Recovery when validation-fold prediction fails")

run_varimpact = function(data, Y, verbose = FALSE, ...) {
  varimpact(Y = Y, data = data, V = 2L,
            Q.library = c("SL.mean", "SL.glm"),
            g.library = c("SL.mean", "SL.glm"),
            verbose = verbose, ...)
}

set.seed(5, "L'Ecuyer-CMRG")
n = 250
X_num = data.frame(x1 = rnorm(n), x2 = rnorm(n))
X_fac = data.frame(f1 = as.factor(sample(c("a", "b", "c"), n, replace = TRUE)),
                   f2 = as.factor(sample(c("a", "b", "c"), n, replace = TRUE)))
Y_bin = rbinom(n, 1, plogis(1.5 * X_num$x1))

future::plan("sequential")

test_that("numerics survive a failed low-level prediction", {
  # Odd calls are the low (control) level of each fold.
  vim = with_failing("apply_tmle_to_validation", function(i) i %% 2 == 1,
                     run_varimpact(X_num, Y_bin))
  expect_s3_class(vim, "varimpact")
})

test_that("numerics survive a failed high-level prediction", {
  vim = with_failing("apply_tmle_to_validation", function(i) i %% 2 == 0,
                     run_varimpact(X_num, Y_bin))
  expect_s3_class(vim, "varimpact")
})

test_that("factors survive failed predictions at either level", {
  vim = with_failing("apply_tmle_to_validation", function(i) i %% 2 == 1,
                     run_varimpact(X_fac, Y_bin))
  expect_s3_class(vim, "varimpact")
  vim = with_failing("apply_tmle_to_validation", function(i) i %% 2 == 0,
                     run_varimpact(X_fac, Y_bin))
  expect_s3_class(vim, "varimpact")
})

test_that("nothing is reported when every validation prediction fails", {
  expect_warning(
    vim <- with_failing("apply_tmle_to_validation", function(i) TRUE,
                        run_varimpact(X_num, Y_bin)),
    "No VIMs could be calculated")
  expect_s3_class(vim, "varimpact")
  expect_true(is.null(vim$results_all) || nrow(vim$results_all) == 0)

  # Factors take a separate path to the same outcome, and when verbose each
  # bin reports that no fold produced rows for it.
  out = capture.output(suppressWarnings(
    vim <- with_failing("apply_tmle_to_validation", function(i) TRUE,
                        run_varimpact(X_fac, Y_bin, verbose = TRUE))))
  expect_s3_class(vim, "varimpact")
  expect_true(any(grepl("no rows are available", out, fixed = TRUE)))
})

test_that("apply_tmle_to_validation is restored afterwards", {
  vim = run_varimpact(X_num, Y_bin)
  expect_true(nrow(vim$results_all) > 0)
})

test_that("the recovery paths report themselves when verbose", {
  for (X in list(X_num, X_fac)) {
    out = capture.output(
      vim <- with_failing("apply_tmle_to_validation", function(i) i %% 2 == 1,
                          run_varimpact(X, Y_bin, verbose = TRUE)))
    expect_s3_class(vim, "varimpact")
    expect_true(any(grepl("prediction on validation failed", out, fixed = TRUE)))
    out = capture.output(
      vim <- with_failing("apply_tmle_to_validation", function(i) i %% 2 == 0,
                          run_varimpact(X, Y_bin, verbose = TRUE)))
    expect_true(any(grepl("prediction on validation failed", out, fixed = TRUE)))
  }
})
