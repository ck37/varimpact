# Regression tests for GitHub issue #32: varimpact() failed whenever Y (or A)
# contained missing values, because observations missing an outcome were only
# ever dropped from the training fold - and then only when fewer than 10 of
# them were missing - so the validation fold always received NA outcomes.
library(varimpact)
library(testthat)

context("Missing outcome values")

future::plan("sequential")

# Simulation from the issue report.
sim_data = function(num_missing_y, N = 200L) {
  set.seed(1, "L'Ecuyer-CMRG")
  num_normal = 4L
  X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
  Y = rbinom(N, 1, plogis(0.2 * X[, 1] + 0.1 * X[, 2] - 0.2 * X[, 3] +
                          0.1 * X[, 3] * X[, 4] - 0.2 * abs(X[, 4])))
  # Add some missing data to X so that imputation is also exercised.
  for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA
  missing_rows = c(4, 6, 7, 8, 11, 15, 20, 21, 28, 32, 72)
  Y[missing_rows[seq_len(num_missing_y)]] = NA
  list(X = X, Y = Y)
}

test_that("varimpact() estimates VIMs when Y has missing values", {
  data = sim_data(11L)
  expect_equal(sum(is.na(data$Y)), 11L)

  vim = varimpact(Y = data$Y, data = data$X)

  expect_s3_class(vim, "varimpact")
  expect_true(nrow(vim$results_all) > 0L)
  expect_false(any(is.na(vim$results_all$Estimate)))
  # A missing outcome must not turn the standard errors into NAs.
  expect_false(any(is.na(vim$results_all[["P-value"]])))
})

test_that("missing outcomes are handled the same above and below 10", {
  # The old code dropped training observations only when fewer than 10 were
  # missing, so 9 and 11 missing outcomes took different code paths.
  for (num_missing in c(9L, 11L)) {
    data = sim_data(num_missing)
    vim = varimpact(Y = data$Y, data = data$X)
    expect_true(nrow(vim$results_all) > 0L,
                info = paste("num_missing =", num_missing))
    expect_false(any(is.na(vim$results_all$Estimate)),
                 info = paste("num_missing =", num_missing))
  }
})

test_that("estimate_tmle2() handles a missing outcome via the delta mechanism", {
  set.seed(2, "L'Ecuyer-CMRG")
  N = 200L
  W = data.frame(w1 = rnorm(N), w2 = rnorm(N))
  A = rbinom(N, 1, plogis(0.5 * W$w1))
  Y = rbinom(N, 1, plogis(0.4 * A + 0.3 * W$w1 - 0.2 * W$w2))
  Y[1:20] = NA
  delta = as.numeric(!is.na(Y))

  result = estimate_tmle2(Y, A, W, family = "binomial", delta = delta,
                          Q.lib = c("SL.mean", "SL.glm"),
                          g.lib = c("SL.mean", "SL.glm"))

  expect_true(is.finite(result$theta))
  # Observations are no longer dropped, so the influence curve covers them all.
  expect_equal(length(result$IC), N)
  expect_false(any(is.na(result$IC)))
  # The missingness model is needed to apply this fit to a validation fold.
  expect_false(is.null(result$g_delta_model))
})

test_that("estimate_tmle2() returns no missingness model without missingness", {
  set.seed(3, "L'Ecuyer-CMRG")
  N = 150L
  W = data.frame(w1 = rnorm(N), w2 = rnorm(N))
  A = rbinom(N, 1, plogis(0.5 * W$w1))
  Y = rbinom(N, 1, plogis(0.4 * A + 0.3 * W$w1))

  result = estimate_tmle2(Y, A, W, family = "binomial",
                          Q.lib = c("SL.mean", "SL.glm"),
                          g.lib = c("SL.mean", "SL.glm"))

  expect_null(result$g_delta_model)
  expect_true(is.finite(result$theta))
})

test_that("validation predictions zero out observations with a missing outcome", {
  set.seed(4, "L'Ecuyer-CMRG")
  N = 200L
  W = data.frame(w1 = rnorm(N), w2 = rnorm(N))
  A = rbinom(N, 1, plogis(0.5 * W$w1))
  Y = rbinom(N, 1, plogis(0.4 * A + 0.3 * W$w1))
  Y[1:20] = NA
  delta = as.numeric(!is.na(Y))

  train = 1:100
  valid = 101:200
  fit = estimate_tmle2(Y[train], A[train], W[train, ], family = "binomial",
                       delta = delta[train],
                       Qbounds = c(0, 1),
                       Q.lib = c("SL.mean", "SL.glm"),
                       g.lib = c("SL.mean", "SL.glm"))

  preds = varimpact:::apply_tmle_to_validation(
    Y = Y[valid], A = A[valid], W = W[valid, ], family = "binomial",
    delta = delta[valid], tmle = fit)

  expect_equal(nrow(preds), length(valid))
  # Nothing may be NA, otherwise the pooled fluctuation and influence curve
  # collapse to NA.
  expect_false(any(is.na(preds)))
  # Missing outcomes get zero weight in the fluctuation.
  expect_true(all(preds$HAW[preds$delta == 0] == 0))
  expect_true(all(preds$HAW[preds$delta == 1 & preds$A == 1] > 0))
})

test_that("create_cv_folds() stratifies a binary outcome that has NAs", {
  set.seed(5, "L'Ecuyer-CMRG")
  Y = c(rep(0, 60), rep(1, 60), rep(NA, 20))
  folds = varimpact:::create_cv_folds(V = 2L, Y = Y)

  expect_equal(length(folds), length(Y))
  expect_false(any(is.na(folds)))
  # Each of the three strata (Y = 0, Y = 1, Y missing) is split evenly.
  counts = table(ifelse(is.na(Y), "missing", Y), folds)
  expect_true(all(counts == c(30, 30, 10)))
})
