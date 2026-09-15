# Issue #8: variable importance estimates for a continuous outcome were reported
# on the internal [0, 1] scale rather than the scale of Y itself.
#
# varimpact() maps a continuous outcome into [0, 1] with
# Y_star = (Y - Qbounds[1]) / diff(Qbounds) before running the CV-TMLE, where
# Qbounds is the observed range of Y widened by 10%. estimate_pooled_results()
# now applies the inverse of that map, so both the treatment-specific means and
# their influence curves come back on the outcome's own scale.

library(varimpact)
library(testthat)

context("Continuous outcome scale transformation")

# Build a minimal fold_results structure of the shape estimate_pooled_results()
# consumes: one element per fold, each holding a $val_preds data frame with the
# columns apply_tmle_to_validation() returns.
make_fold_results <- function(n_per_fold = 60, n_folds = 2, seed = 7) {
  set.seed(seed)
  lapply(seq_len(n_folds), function(fold) {
    A <- rbinom(n_per_fold, 1, 0.5)
    g1W_hat <- runif(n_per_fold, 0.3, 0.7)
    Q_hat <- runif(n_per_fold, 0.2, 0.8)
    # Y_star is already on the [0, 1] scale here, as it is in real fold results.
    Y_star <- runif(n_per_fold)
    gDelta_hat <- rep(1, n_per_fold)
    gAW_total <- g1W_hat * gDelta_hat
    H1W <- 1 / gAW_total
    list(val_preds = data.frame(Y_star = Y_star,
                                A = A,
                                Q_hat = Q_hat,
                                g1W_hat = g1W_hat,
                                gDelta_hat = gDelta_hat,
                                gAW_total = gAW_total,
                                delta = rep(1, n_per_fold),
                                H1W = H1W,
                                HAW = A * H1W))
  })
}

test_that("estimate_pooled_results() rescales thetas and ICs by exactly Qbounds", {
  fold_results <- make_fold_results()

  # The internal [0, 1] scale: Qbounds = c(0, 1) must be the identity.
  unit <- varimpact:::estimate_pooled_results(fold_results, verbose = FALSE, Qbounds = c(0, 1))
  default <- varimpact:::estimate_pooled_results(fold_results, verbose = FALSE)
  expect_equal(unit$thetas, default$thetas)

  # An arbitrary continuous-outcome range.
  Qbounds <- c(-12.5, 37.5)
  scaled <- varimpact:::estimate_pooled_results(fold_results, verbose = FALSE, Qbounds = Qbounds)

  # Rescaling happens after the fluctuation, so epsilon is untouched.
  expect_equal(scaled$epsilon, unit$epsilon)

  # theta_original = theta_star * diff(Qbounds) + Qbounds[1]
  expect_equal(as.vector(scaled$thetas),
               as.vector(unit$thetas) * diff(Qbounds) + Qbounds[1])

  # IC_original = diff(Qbounds) * IC_star: the location shift cancels, because
  # both terms of the influence curve are differences.
  for (fold in seq_along(unit$influence_curves)) {
    expect_equal(as.vector(scaled$influence_curves[[fold]]),
                 as.vector(unit$influence_curves[[fold]]) * diff(Qbounds))
  }
})

test_that("estimate_pooled_results() ignores a malformed Qbounds", {
  fold_results <- make_fold_results()
  unit <- varimpact:::estimate_pooled_results(fold_results, verbose = FALSE, Qbounds = c(0, 1))

  expect_equal(varimpact:::estimate_pooled_results(fold_results, Qbounds = NULL)$thetas,
               unit$thetas)
  expect_equal(varimpact:::estimate_pooled_results(fold_results, Qbounds = 4)$thetas,
               unit$thetas)
})

test_that("continuous outcomes produce results on the original scale", {
  # Before the fix to the double plogis() in estimate_tmle2(), every fold of a
  # gaussian run was discarded ("min and max level are the same") and
  # results_all came back NULL, so this is also a regression test for that.
  set.seed(42, "L'Ecuyer-CMRG")
  future::plan("sequential")

  N <- 300
  X <- as.data.frame(matrix(rnorm(N * 2), N, 2))
  colnames(X) <- paste0("X", 1:2)

  # True ATE contrast on the original scale is roughly 4.5 for X1 and 3.0 for X2
  # per unit of X, i.e. far outside [0, 1].
  Y <- 50 + 15 * (0.3 * X[, 1] + 0.2 * X[, 2]) + rnorm(N, 0, 3)
  expect_gt(min(Y), 10)

  vim <- suppressWarnings(
    varimpact(Y = Y, data = X, family = "gaussian", V = 2L, verbose = FALSE,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"),
              bins_numeric = 3L))

  expect_s3_class(vim, "varimpact")
  expect_s3_class(vim$results_all, "data.frame")
  expect_equal(nrow(vim$results_all), 2)

  estimates <- vim$results_all$Estimate
  expect_false(any(is.na(estimates)))

  # The headline check for issue #8: estimates are on Y's scale, not [0, 1].
  # A median split of a standard normal separates the bin means by roughly 1.6,
  # so the X1 contrast should land near 15 * 0.3 * 1.6 = 7.2 and the X2 contrast
  # near 15 * 0.2 * 1.6 = 4.8. Both would be under 1 on the old scale.
  expect_true(all(estimates > 1))
  expect_true(all(estimates < diff(range(Y))))
  expect_equal(estimates[rownames(vim$results_all) == "X1"], 7.2, tolerance = 0.4)
  expect_equal(estimates[rownames(vim$results_all) == "X2"], 4.8, tolerance = 0.4)

  # Confidence intervals are on the same scale, so they bracket the estimate.
  ci <- vim$results_all$CI95
  expect_true(all(!is.na(ci)))
})

test_that("binary outcomes are unaffected by the rescaling", {
  # Qbounds is c(0, 1) for a binary outcome, so every transformation added for
  # issue #8 is the identity and estimates stay within [-1, 1].
  set.seed(3, "L'Ecuyer-CMRG")
  future::plan("sequential")

  N <- 200
  X <- as.data.frame(matrix(rnorm(N * 2), N, 2))
  colnames(X) <- c("X1", "X2")
  Y <- rbinom(N, 1, plogis(0.5 * X[, 1] - 0.3 * X[, 2]))

  vim <- suppressWarnings(
    varimpact(Y = Y, data = X, family = "binomial", V = 2L, verbose = FALSE,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm")))

  expect_s3_class(vim, "varimpact")
  expect_s3_class(vim$results_all, "data.frame")

  estimates <- vim$results_all$Estimate
  estimates <- estimates[!is.na(estimates)]
  expect_gt(length(estimates), 0)
  expect_true(all(estimates >= -1 & estimates <= 1))
})
