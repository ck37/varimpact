library(varimpact)
library(testthat)

# estimate_pooled_results() runs the pooled fluctuation step over the
# validation predictions from every fold. varimpact() only ever hands it
# well-formed input, so its guards are tested here directly, with fold
# results built from a real training-fold fit.

context("estimate_pooled_results")

set.seed(9, "L'Ecuyer-CMRG")
n = 160
W = data.frame(w1 = rnorm(n), w2 = rnorm(n))
A = rbinom(n, 1, plogis(0.4 * W$w1))
Y = rbinom(n, 1, plogis(0.5 * W$w1 + A))
lib = c("SL.mean", "SL.glm")

fit = estimate_tmle2(Y = Y, A = A, W = W, family = "binomial",
                     Q.lib = lib, g.lib = lib, V = 2, Qbounds = c(0, 1))

# Two "validation folds": the two halves of the sample.
halves = split(seq_len(n), rep(1:2, length.out = n))
fold_results = lapply(halves, function(idx) {
  list(val_preds = varimpact:::apply_tmle_to_validation(
    Y[idx], A[idx], W[idx, ], "binomial", delta = rep(1, length(idx)), tmle = fit))
})

test_that("pooling well-formed folds yields one estimate per fold", {
  pooled = varimpact:::estimate_pooled_results(fold_results)
  expect_length(pooled$thetas, 2L)
  expect_false(any(is.na(pooled$thetas)))
  expect_true(all(pooled$thetas >= 0 & pooled$thetas <= 1))
  expect_length(pooled$influence_curves, 2L)
  expect_equal(unname(lengths(pooled$influence_curves)), unname(lengths(halves)))
  expect_true(is.numeric(pooled$epsilon))
})

test_that("verbose reports the fluctuation step", {
  expect_output(varimpact:::estimate_pooled_results(fold_results, verbose = TRUE),
                "Estimating epsilon")
})

test_that("a fold that produced nothing leaves its slot empty", {
  one_missing = fold_results
  one_missing[[2]]$val_preds = NULL
  pooled = varimpact:::estimate_pooled_results(one_missing)
  expect_length(pooled$thetas, 2L)
  expect_false(is.na(pooled$thetas[1]))
  expect_true(is.na(pooled$thetas[2]))
  expect_null(pooled$influence_curves[[2]])
})

test_that("when every fold failed the placeholders come back", {
  none = list(list(val_preds = NULL), list(val_preds = data.frame()))
  expect_output(pooled <- varimpact:::estimate_pooled_results(none, verbose = TRUE),
                "every fold failed")
  expect_null(pooled$thetas)
  expect_null(pooled$influence_curves)
  expect_null(pooled$epsilon)
})

test_that("predicted Q_hat outside [0, 1] is an error, not a debugger prompt", {
  # This used to call browser(), which is a no-op when non-interactive and
  # let the function continue with an invalid state.
  bad = fold_results
  bad[[1]]$val_preds$Q_hat[1] = 1.5
  expect_error(varimpact:::estimate_pooled_results(bad),
               "Q_hat values must lie in \\[0, 1\\]")
})

test_that("only the logistic fluctuation is supported", {
  expect_error(varimpact:::estimate_pooled_results(fold_results, fluctuation = "linear"),
               "Only support logistic fluctuation")
})

test_that("val_preds without a delta column are treated as fully observed", {
  no_delta = lapply(fold_results, function(fold) {
    fold$val_preds$delta = NULL
    fold
  })
  pooled = varimpact:::estimate_pooled_results(no_delta)
  expect_equal(pooled$thetas, varimpact:::estimate_pooled_results(fold_results)$thetas)
})
