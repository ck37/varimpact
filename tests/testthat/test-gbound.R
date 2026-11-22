library(varimpact)
library(testthat)

context("gbound: truncation of the estimated propensity score")

future::plan("sequential")

# One numeric and one factor variable, so that vim_numerics() and
# vim_factors() both fit TMLEs. Level "c" is held by about 10% of
# observations.
set.seed(5, "L'Ecuyer-CMRG")
n = 200
X = data.frame(num = rnorm(n),
               fac = factor(sample(c("a", "b", "c"), n, replace = TRUE,
                                   prob = c(0.6, 0.3, 0.1))))
Y = rbinom(n, 1, plogis(-0.3 + 0.8 * X$num + 0.8 * (X$fac == "c")))
lib = c("SL.mean", "SL.glm")

run = function(...) {
  set.seed(2, "L'Ecuyer-CMRG")
  varimpact(Y = Y, data = X, V = 2L, Q.library = lib, g.library = lib, ...)
}

# Run expr, recording the gbound that each estimate_tmle2() call receives,
# with the function that called it, and the bounds that each
# apply_tmle_to_validation() call applies (taken from its training fit).
record_bounds = function(expr) {
  seen = new.env()
  seen$fits = data.frame(caller = character(0), gbound = numeric(0))
  seen$validation = list()
  caller = function() {
    fns = vapply(sys.calls(), function(cl) deparse(cl[[1]])[1], character(1))
    intersect(c("vim_factors", "vim_numerics"), fns)[1]
  }
  value = with_replaced("estimate_tmle2", function(original) {
    function(...) {
      gbound = list(...)[["gbound"]]
      seen$fits = rbind(seen$fits, data.frame(
        caller = caller(), gbound = if (is.null(gbound)) NA_real_ else gbound))
      original(...)
    }
  }, with_replaced("apply_tmle_to_validation", function(original) {
    function(Y, A, W, family, delta = rep(1, length(Y)), tmle, ...) {
      seen$validation[[length(seen$validation) + 1L]] = tmle$gbounds
      original(Y, A, W, family, delta, tmle, ...)
    }
  }, expr))
  list(value = value, fits = seen$fits,
       lower = vapply(seen$validation, `[`, numeric(1), 1L))
}

test_that("gbound defaults to 0.025 for factors and numerics alike", {
  r = record_bounds(run())
  expect_identical(r$value$gbound, 0.025)
  expect_setequal(r$fits$caller, c("vim_factors", "vim_numerics"))
  expect_true(all(r$fits$gbound == 0.025))
  expect_gt(length(r$lower), 0L)
  expect_true(all(r$lower == 0.025))
})

test_that("a chosen gbound reaches every fit and every validation fold", {
  r = record_bounds(run(gbound = 0.01))
  expect_identical(r$value$gbound, 0.01)
  expect_setequal(r$fits$caller, c("vim_factors", "vim_numerics"))
  expect_true(all(r$fits$gbound == 0.01))
  expect_gt(length(r$lower), 0L)
  expect_true(all(r$lower == 0.01))
})

test_that("a bound that binds moves the estimates; the default does not", {
  default = run()
  expect_identical(run(gbound = 0.025)$results_raw, default$results_raw)
  # A 10% level has estimated probabilities below 0.2, which a bound of 0.2
  # raises, so the estimates change.
  expect_false(identical(run(gbound = 0.2)$results_raw, default$results_raw))
})
