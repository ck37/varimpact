# Regression tests for the pooled CV-TMLE fluctuation in estimate_pooled_results().
#
# The formula used to say stats::offset(logit_Q_hat). Namespace-qualifying
# offset() stops R's terms machinery from recognizing it as the formula's
# offset special, so logit_Q_hat was fit as an ordinary covariate: epsilon came
# back with two elements instead of one, and "epsilon * data$H1W" then recycled
# them alternately across observations.
#
# These use a fixed synthetic pair of validation folds rather than a varimpact()
# run, so they are deterministic and fast. An end-to-end check on null data
# separates the two code paths only in distribution, not per seed, which would
# make for a flaky test.
library(varimpact)
library(testthat)

context("Pooled fluctuation")

# Build one fold in the shape apply_tmle_to_validation() returns.
make_fold = function(n, seed) {
  set.seed(seed, "L'Ecuyer-CMRG")
  g1W_hat = runif(n, 0.2, 0.8)
  gDelta_hat = runif(n, 0.8, 1.0)
  gAW_total = g1W_hat * gDelta_hat
  A = rbinom(n, 1, g1W_hat)
  delta = rbinom(n, 1, gDelta_hat)
  Q_hat = runif(n, 0.2, 0.8)
  # Q is well specified: Y is drawn from Q_hat itself.
  Y_star = rbinom(n, 1, Q_hat)
  # Missing outcomes carry zero weight, so their value is a placeholder.
  Y_star[delta == 0] = 0
  H1W = 1 / gAW_total

  list(val_preds = data.frame(Y_star = Y_star, A = A, Q_hat = Q_hat,
                              g1W_hat = g1W_hat, gDelta_hat = gDelta_hat,
                              gAW_total = gAW_total, delta = delta,
                              H1W = H1W, HAW = A * delta * H1W))
}

# Odd row counts on purpose: an even total hid the recycling behind a silent
# multiple, so only odd counts produced the "longer object length" warning.
pooled_fit = function() {
  folds = list(make_fold(99L, 1L), make_fold(101L, 2L))
  list(result = varimpact:::estimate_pooled_results(folds), folds = folds)
}

test_that("the fluctuation estimates exactly one coefficient", {
  epsilon = pooled_fit()$result$epsilon

  # logit_Q_hat is an offset, not a covariate, so it must not get a
  # coefficient of its own. Two elements here means offset() was not honored
  # and epsilon will recycle across observations downstream.
  expect_length(epsilon, 1L)
  expect_identical(names(epsilon), "HAW")
  expect_true(is.finite(epsilon))
})

test_that("per-fold influence curves are finite and centered near zero", {
  fit = pooled_fit()
  curves = fit$result$influence_curves

  expect_length(curves, 2L)
  expect_equal(unname(vapply(curves, length, integer(1))), c(99L, 101L))
  expect_true(all(vapply(curves, function(x) all(is.finite(x)), logical(1))))

  # A correctly targeted influence curve has an empirical mean near zero, up
  # to Monte Carlo error of order 1/sqrt(n) (about 0.1 at n = 100). The
  # recycled epsilon pushed this to 0.35.
  ic_means = vapply(curves, mean, numeric(1))
  expect_true(all(abs(ic_means) < 0.15))
})

test_that("targeting does not move the estimate far from initial Q", {
  fit = pooled_fit()

  # Q is well specified here, so the fluctuation should be small and each
  # fold's estimate should land near the mean of the initial Q_hat.
  mean_q = mean(fit$folds[[1]]$val_preds$Q_hat)
  expect_true(all(abs(fit$result$thetas - mean_q) < 0.15))
})
