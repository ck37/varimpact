# The weighted (intercept-only) fluctuation, option varimpact.weighted_fluctuation.
#
# The covariate form updates each unit by epsilon * H1W. Under a near-violation
# of positivity many units sit at the g bound (H1W = 1 / 0.025 = 40), and a
# moderate epsilon fitted on the treated units sends their counterfactual Q* to 1.
# The weighted form fits epsilon as the intercept of a submodel weighted by the
# clever covariate, so every unit moves by the same epsilon. Both solve the same
# efficient influence curve equation.
library(varimpact)
library(testthat)

context("Weighted fluctuation")

# One validation fold in the shape apply_tmle_to_validation() returns. With
# near_violation = TRUE, a third of the units have a treatment probability at
# the 0.025 bound and are untreated, while the treated units have a moderate
# propensity and an event rate five times the initial Q.
make_fold = function(n, seed, near_violation = FALSE) {
  set.seed(seed, "L'Ecuyer-CMRG")
  g1W_hat = if (near_violation) ifelse(seq_len(n) %% 3 == 0, 0.025, 0.5) else runif(n, 0.2, 0.8)
  gDelta_hat = rep(1, n)
  gAW_total = g1W_hat * gDelta_hat
  A = if (near_violation) as.integer(g1W_hat > 0.025 & seq_len(n) %% 2 == 0) else rbinom(n, 1, g1W_hat)
  delta = rep(1L, n)
  Q_hat = if (near_violation) rep(0.01, n) else runif(n, 0.2, 0.8)
  Y_star = if (near_violation) as.numeric(A == 1 & seq_len(n) %% 20 == 0) else rbinom(n, 1, Q_hat)
  H1W = 1 / gAW_total
  list(val_preds = data.frame(Y_star = Y_star, A = A, Q_hat = Q_hat, g1W_hat = g1W_hat,
                              gDelta_hat = gDelta_hat, gAW_total = gAW_total, delta = delta,
                              H1W = H1W, HAW = A * delta * H1W))
}

pooled = function(weighted, near_violation = FALSE) {
  folds = list(make_fold(600L, 1L, near_violation), make_fold(600L, 2L, near_violation))
  varimpact:::estimate_pooled_results(folds, weighted_fluctuation = weighted)
}

test_that("the weighted fluctuation solves the influence curve equation", {
  folds = list(make_fold(99L, 1L), make_fold(101L, 2L))
  res = varimpact:::estimate_pooled_results(folds, weighted_fluctuation = TRUE)
  expect_length(res$epsilon, 1L)
  d = do.call(rbind, lapply(folds, `[[`, "val_preds"))
  q_star = plogis(qlogis(d$Q_hat) + res$epsilon)
  # sum over observed treated units of H1W * (Y - Q*) is zero at the fitted epsilon.
  expect_lt(abs(sum(d$HAW * (d$Y_star - q_star))), 1e-6)
})

test_that("on well-behaved data the two forms agree", {
  w = pooled(TRUE); c = pooled(FALSE)
  expect_equal(unname(w$thetas), unname(c$thetas), tolerance = 0.02)
})

test_that("under a near-violation the covariate form saturates and the weighted form does not", {
  c = pooled(FALSE, near_violation = TRUE)
  w = pooled(TRUE, near_violation = TRUE)
  # Treated event rate is 0.1 against an initial Q of 0.01; a third of the units
  # sit at the bound with H1W = 40.
  expect_gt(min(c$thetas), 0.3)              # Q* driven toward 1 for the units at the bound
  expect_lt(max(w$thetas), 0.15)             # one shift for everyone: near the treated event rate
  expect_true(all(is.finite(w$influence_curves[[1]])))
})

test_that("the weighted form is the default, and the option switches it", {
  folds = list(make_fold(99L, 1L), make_fold(101L, 2L))
  old = options(varimpact.weighted_fluctuation = NULL); on.exit(options(old))
  expect_identical(varimpact:::estimate_pooled_results(folds)$epsilon,
                   varimpact:::estimate_pooled_results(folds, weighted_fluctuation = TRUE)$epsilon)
  options(varimpact.weighted_fluctuation = FALSE)
  expect_identical(varimpact:::estimate_pooled_results(folds)$epsilon,
                   varimpact:::estimate_pooled_results(folds, weighted_fluctuation = FALSE)$epsilon)
})

# The training-fold TMLE in estimate_tmle2(): with a well-specified glm library
# both forms solve the influence curve equation for E[Y | A = 1, W], so the
# returned influence curve has mean zero and the two estimates agree.
sim_point = function(n = 2000L, seed = 11L) {
  set.seed(seed)
  W = data.frame(w1 = rnorm(n), w2 = rbinom(n, 1, 0.4))
  A = rbinom(n, 1, plogis(-0.3 + 0.6 * W$w1))
  Y = rbinom(n, 1, plogis(-1.5 + 0.7 * A + 0.5 * W$w1 - 0.4 * W$w2))
  list(Y = Y, A = A, W = W)
}

test_that("estimate_tmle2's weighted fluctuation solves the influence curve equation", {
  s = sim_point()
  w = varimpact:::estimate_tmle2(s$Y, s$A, s$W, family = "binomial", Q.lib = "SL.glm", g.lib = "SL.glm",
                                 V = 2, weighted_fluctuation = TRUE)
  c = varimpact:::estimate_tmle2(s$Y, s$A, s$W, family = "binomial", Q.lib = "SL.glm", g.lib = "SL.glm",
                                 V = 2, weighted_fluctuation = FALSE)
  expect_lt(abs(mean(w$IC)), 1e-6)
  expect_true(w$theta > 0 && w$theta < 1)
  expect_equal(w$theta, c$theta, tolerance = 0.01)
})

test_that("varimpact() takes the argument, records it, and the two forms agree on well-behaved data", {
  s = sim_point(n = 1500L, seed = 12L)
  X = data.frame(A = s$A, s$W)
  run = function(...) {
    set.seed(5)
    suppressMessages(varimpact(Y = s$Y, data = X, A_names = "A", V = 2L, Q.library = c("SL.mean", "SL.glm"),
                               g.library = c("SL.mean", "SL.glm"), family = "binomial", bins_numeric = 2L,
                               quantile_probs_factor = NULL, quantile_probs_numeric = NULL, verbose = FALSE,
                               parallel = FALSE, ...))
  }
  old = options(varimpact.weighted_fluctuation = NULL); on.exit(options(old))
  default = run()
  covariate = run(weighted_fluctuation = FALSE)
  expect_true(default$weighted_fluctuation)
  expect_false(covariate$weighted_fluctuation)
  expect_true(is.finite(default$results_all$Estimate))
  expect_equal(default$results_all$Estimate, covariate$results_all$Estimate, tolerance = 0.03)
  # The argument overrides the option.
  options(varimpact.weighted_fluctuation = FALSE)
  expect_true(run(weighted_fluctuation = TRUE)$weighted_fluctuation)
})

test_that("varimpact() rejects a weighted_fluctuation that is not TRUE or FALSE", {
  s = sim_point(n = 200L, seed = 13L)
  X = data.frame(A = s$A, s$W)
  for (bad in list(NA, "yes", c(TRUE, FALSE)))
    expect_error(varimpact(Y = s$Y, data = X, weighted_fluctuation = bad), "TRUE or FALSE")
})
