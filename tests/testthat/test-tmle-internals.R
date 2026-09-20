library(varimpact)
library(SuperLearner)
library(testthat)

# Direct tests of the TMLE building blocks that varimpact() only ever calls
# one way. Each has argument checks and fallbacks that the end-to-end tests
# never reach: tmle_init_stage1() (the outcome transform), tmle_estimate_q()
# (the initial Q fit), estimate_tmle2() (the training-fold fit) and
# apply_tmle_to_validation() (the validation-fold predictions).

context("TMLE internals")

set.seed(7, "L'Ecuyer-CMRG")
n = 120
W = data.frame(w1 = rnorm(n), w2 = rnorm(n))
A = rbinom(n, 1, plogis(0.4 * W$w1))
Y_cont = 1 + W$w1 - 0.5 * W$w2 + A + rnorm(n)
Y_bin = rbinom(n, 1, plogis(0.5 * W$w1 + A))
Delta = rep(1, n)
lib = c("SL.mean", "SL.glm")

# ---- tmle_init_stage1 -------------------------------------------------------

test_that("a continuous outcome is mapped onto [0, 1] with widened bounds", {
  res = varimpact:::tmle_init_stage1(Y = Y_cont, A = A, Q = NULL, Delta = Delta,
                                     Qbounds = NULL, alpha = 0.995,
                                     maptoYstar = TRUE, family = "gaussian")
  expect_equal(res$ab, range(Y_cont))
  expect_equal(range(res$Ystar), c(0, 1))
  expect_equal(res$Ystar, (Y_cont - min(Y_cont)) / diff(range(Y_cont)))
  # The returned Qbounds are the fluctuation bounds, alpha and 1 - alpha.
  expect_equal(sort(res$Qbounds), c(0.005, 0.995))
})

test_that("without the mapping the outcome and bounds are left alone", {
  res = varimpact:::tmle_init_stage1(Y = Y_cont, A = A, Q = NULL, Delta = Delta,
                                     Qbounds = NULL, alpha = 0.995,
                                     maptoYstar = FALSE, family = "gaussian")
  expect_identical(res$Ystar, Y_cont)
  expect_equal(res$Qbounds, c(-Inf, Inf))
  expect_equal(res$ab, c(0, 1))
})

test_that("an alpha outside (0, 1) is reset with a warning", {
  expect_warning(
    res <- varimpact:::tmle_init_stage1(Y = Y_cont, A = A, Q = NULL, Delta = Delta,
                                        Qbounds = NULL, alpha = 1.5,
                                        maptoYstar = TRUE, family = "gaussian"),
    "alpha must be between 0 and 1")
  expect_equal(sort(res$Qbounds), c(0.005, 0.995))
})

test_that("user-supplied Q values are expanded to QAW, Q0W, Q1W and rescaled", {
  Q = cbind(Q0W = plogis(W$w1), Q1W = plogis(W$w1 + 1))
  res = varimpact:::tmle_init_stage1(Y = Y_bin, A = A, Q = Q, Delta = Delta,
                                     Qbounds = NULL, alpha = 0.995,
                                     maptoYstar = TRUE, family = "binomial")
  expect_equal(colnames(res$Q), c("QAW", "Q0W", "Q1W"))
  # Binary outcome: bounds are [0, 1], so the rescaling is the identity and
  # QAW is the row's own-arm prediction.
  expect_equal(res$Q[, "QAW"], (1 - A) * Q[, "Q0W"] + A * Q[, "Q1W"])
  expect_equal(res$ab, c(0, 1))
})

# ---- tmle_estimate_q --------------------------------------------------------

test_that("Qbounds must be supplied", {
  expect_error(
    tmle_estimate_q(Y = Y_bin, A = A, W = W, Delta = Delta, Qbounds = NULL,
                    maptoYstar = TRUE, SL.library = lib, family = "binomial"),
    "Qbounds must be defined")
})

test_that("a Qform of Y ~ . fits a main-terms glm instead of SuperLearner", {
  q = tmle_estimate_q(Y = Y_bin, A = A, W = W, Delta = Delta, Qbounds = c(0, 1),
                      Qform = "Y ~ .", maptoYstar = FALSE, SL.library = lib,
                      family = "binomial")
  expect_equal(q$type, "glm, user-supplied model")
  expect_s3_class(q$model, "glm")
  expect_equal(colnames(q$Q), c("QAW", "Q0W", "Q1W"))
  # Binomial Q is returned on the logit scale; undoing it recovers the glm's
  # own predictions under A = 1.
  expected = predict(q$model, newdata = data.frame(Y = Y_bin, Z = 0, A = 1, W),
                     type = "response")
  expect_equal(unname(plogis(q$Q[, "Q1W"])), unname(expected))
  expect_equal(q$family, "binomial")
})

test_that("a poisson family returns Q on the log scale", {
  Y_count = rpois(n, exp(0.2 * W$w1 + 0.5 * A))
  q = tmle_estimate_q(Y = Y_count, A = A, W = W, Delta = Delta,
                      Qbounds = c(0, Inf), Qform = "Y ~ .", maptoYstar = FALSE,
                      SL.library = lib, family = "poisson")
  expect_equal(q$family, "poisson")
  expected = predict(q$model, newdata = data.frame(Y = Y_count, Z = 0, A = 0, W),
                     type = "response")
  expect_equal(unname(exp(q$Q[, "Q0W"])), unname(expected))
})

test_that("the SuperLearner path reports itself and keeps the library", {
  expect_output(
    q <- tmle_estimate_q(Y = Y_bin, A = A, W = W, Delta = Delta, Qbounds = c(0, 1),
                         maptoYstar = TRUE, SL.library = lib, family = "binomial",
                         V = 2, verbose = TRUE),
    "using SuperLearner")
  expect_equal(q$type, "SuperLearner")
  expect_equal(q$SL.library, lib)
  expect_s3_class(q$model, "SuperLearner")
})

test_that("a SuperLearner failure is an error, not a silent fallback", {
  expect_error(
    tmle_estimate_q(Y = Y_bin, A = A, W = W, Delta = Delta, Qbounds = c(0, 1),
                    maptoYstar = TRUE, SL.library = "SL.no_such_learner",
                    family = "binomial", V = 2),
    "Super Learner failed")
})

# ---- estimate_tmle2 ---------------------------------------------------------

test_that("estimate_tmle2 validates its family and the shape of W", {
  expect_error(estimate_tmle2(Y = Y_bin, A = A, W = W, family = "poisson",
                              Q.lib = lib, g.lib = lib),
               "family must be either")
  expect_error(estimate_tmle2(Y = Y_bin, A = A, W = W$w1, family = "binomial",
                              Q.lib = lib, g.lib = lib),
               "W should have two dimensions")
})

test_that("an explicit NULL delta means every observation is observed", {
  fit = estimate_tmle2(Y = Y_bin, A = A, W = W, family = "binomial",
                       delta = NULL, Q.lib = lib, g.lib = lib, V = 2)
  expect_true(is.numeric(fit$theta))
  expect_length(fit$IC, n)
  expect_null(fit$g_delta_model)
})

test_that("missing values in W are reported before the fit is attempted", {
  W_na = W
  W_na$w1[1:3] = NA
  # SuperLearner refuses missing predictors, so the fit itself fails; what is
  # being checked is that the report of where the NAs are comes first.
  expect_output(
    try(estimate_tmle2(Y = Y_bin, A = A, W = W_na, family = "binomial",
                       Q.lib = lib, g.lib = lib, V = 2), silent = TRUE),
    "found 3 NAs in W")
})

# ---- apply_tmle_to_validation -----------------------------------------------

# varimpact() passes Qbounds down explicitly for a continuous outcome: the
# range of Y widened by 10% at each end. apply_tmle_to_validation() reads them
# back off the fit, so the fit has to carry them.
Qbounds = range(Y_cont) + 0.1 * c(-abs(min(Y_cont)), abs(max(Y_cont)))
fit = estimate_tmle2(Y = Y_cont, A = A, W = W, family = "gaussian",
                     Q.lib = lib, g.lib = lib, V = 2, Qbounds = Qbounds)

test_that("validation predictions come back one row per observation", {
  val = varimpact:::apply_tmle_to_validation(Y_cont, A, W, "gaussian",
                                             delta = Delta, tmle = fit)
  expect_s3_class(val, "data.frame")
  expect_equal(nrow(val), n)
  expect_true(all(c("Y_star", "A", "Q_hat", "g1W_hat", "HAW", "delta") %in% names(val)))
  # The clever covariate is zero off treatment and Q_hat respects its bounds.
  expect_true(all(val$HAW[val$A == 0] == 0))
  expect_true(all(val$Q_hat >= 0 & val$Q_hat <= 1))
  expect_true(all(val$Y_star >= 0 & val$Y_star <= 1))
})

test_that("observations flagged missing get zero weight and a placeholder outcome", {
  Y_miss = Y_cont
  Y_miss[1:5] = NA
  delta = as.numeric(!is.na(Y_miss))
  val = varimpact:::apply_tmle_to_validation(Y_miss, A, W, "gaussian",
                                             delta = delta, tmle = fit)
  expect_equal(val$HAW[1:5], rep(0, 5))
  expect_equal(val$Y_star[1:5], rep(0, 5))
  expect_false(any(is.na(val$Y_star)))
})

test_that("a missing outcome that is flagged as observed is refused", {
  Y_miss = Y_cont
  Y_miss[1] = NA
  expect_error(
    varimpact:::apply_tmle_to_validation(Y_miss, A, W, "gaussian",
                                         delta = Delta, tmle = fit),
    "delta must be 0")
})

test_that("a validation set with nothing observed is refused", {
  expect_error(
    varimpact:::apply_tmle_to_validation(Y_cont, A, W, "gaussian",
                                         delta = rep(0, n), tmle = fit),
    "no observations have both Y and A observed")
})

test_that("outcomes outside the training range are refused rather than clipped", {
  # Qbounds come from the training fold's range of Y (widened by 10%), so a
  # validation outcome far beyond it cannot be mapped into [0, 1].
  expect_error(
    varimpact:::apply_tmle_to_validation(Y_cont * 100, A, W, "gaussian",
                                         delta = Delta, tmle = fit),
    "0 <= y <= 1")
})

test_that("verbose reports the outcome mapping", {
  expect_output(
    varimpact:::apply_tmle_to_validation(Y_cont, A, W, "gaussian",
                                         delta = Delta, tmle = fit, verbose = TRUE),
    "Mapped Y to Y_star")
})

# ---- the paths varimpact() never takes --------------------------------------

test_that("estimate_tmle2 reports missing values in A and Y before fitting", {
  A_na = A
  A_na[1:2] = NA
  expect_output(
    try(estimate_tmle2(Y = Y_bin, A = A_na, W = W, family = "binomial",
                       Q.lib = lib, g.lib = lib, V = 2), silent = TRUE),
    "found 2 NAs in A")
  Y_na = Y_bin
  Y_na[1:4] = NA
  expect_output(
    try(estimate_tmle2(Y = Y_na, A = A, W = W, family = "binomial",
                       Q.lib = lib, g.lib = lib, V = 2), silent = TRUE),
    "found 4 NAs in Y")
})

test_that("estimate_tmle2 reduces the CV folds when A or delta has a small cell", {
  # Three treated observations and ten folds: the propensity model cannot be
  # cross-validated ten ways, so the folds are reduced to the cell size.
  A_rare = c(rep(1, 3), rep(0, n - 3))
  expect_output(
    estimate_tmle2(Y = Y_cont, A = A_rare, W = W, family = "gaussian",
                   Q.lib = "SL.mean", g.lib = "SL.mean", V = 10, verbose = TRUE),
    "Reducing folds to 3")
  # The same for the missingness mechanism.
  delta_rare = c(rep(0, 3), rep(1, n - 3))
  expect_output(
    estimate_tmle2(Y = Y_cont, A = A, W = W, family = "gaussian", delta = delta_rare,
                   Q.lib = "SL.mean", g.lib = "SL.mean", V = 10, verbose = TRUE),
    "Delta's minimum cell size 3")
})

test_that("tmle_estimate_q falls back to a main-terms glm when supplied Q values are unusable", {
  Q_bad = matrix(NA_real_, nrow = n, ncol = 3,
                 dimnames = list(NULL, c("QAW", "Q0W", "Q1W")))
  expect_output(
    q <- tmle_estimate_q(Y = Y_bin, A = A, W = W, Delta = Delta, Q = Q_bad,
                         Qbounds = c(0, 1), maptoYstar = FALSE, SL.library = lib,
                         family = "binomial", verbose = TRUE),
    "main terms regression")
  expect_equal(q$type, "glm, main terms model")
  expect_s3_class(q$model, "glm")
  expect_false(any(is.na(q$Q)))
})

test_that("apply_tmle_to_validation reports which prediction failed", {
  broken_q = fit
  broken_q$q_model = list(not = "a model")
  expect_error(
    suppressWarnings(varimpact:::apply_tmle_to_validation(
      Y_cont, A, W, "gaussian", delta = Delta, tmle = broken_q)),
    "failed during prediction of Q\\(1, W\\)")

  broken_delta = fit
  broken_delta$g_delta_model = list(not = "a model")
  expect_error(
    varimpact:::apply_tmle_to_validation(Y_cont, A, W, "gaussian",
                                         delta = Delta, tmle = broken_delta),
    "failed during prediction of g.Delta")
})
