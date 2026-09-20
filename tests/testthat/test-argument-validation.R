library(varimpact)
library(testthat)

context("Argument validation in varimpact()")

fit = function(...) {
  set.seed(1, "L'Ecuyer-CMRG")
  n = 150
  X = data.frame(a = rnorm(n), b = rnorm(n))
  X$a[3] = NA
  Y = rbinom(n, 1, 0.5)
  varimpact(Y = Y, data = X, V = 2L,
            Q.library = c("SL.mean", "SL.glm"),
            g.library = c("SL.mean", "SL.glm"),
            verbose = FALSE, ...)
}

future::plan("sequential")

test_that("an unrecognized impute method is reported as such", {
  # Previously this fell through every branch in process_numerics() and
  # surfaced as "sum(is.na(data.numW)) == 0 is not TRUE", which named neither
  # the argument nor the valid values.
  expect_error(fit(impute = "nonesuch"), "should be one of")
  expect_error(fit(impute = "nonesuch"), "median")
})

test_that("mean imputation still reports that it is not implemented", {
  # "mean" is a recognized value, so it reaches its own message rather than
  # being reported as an invalid one.
  expect_error(fit(impute = "mean"), "not implemented")
})

test_that("the supported impute methods are accepted", {
  for (method in c("median", "knn", "zero")) {
    expect_s3_class(fit(impute = method), "varimpact")
  }
})

test_that("a binomial outcome must lie in [0, 1]", {
  X = data.frame(a = rnorm(30), b = rnorm(30))
  expect_error(varimpact(Y = rpois(30, 3), data = X, family = "binomial"),
               "bounded by \\[0, 1\\]")
})

test_that("only binomial and gaussian families are accepted", {
  X = data.frame(a = rnorm(30), b = rnorm(30))
  expect_error(varimpact(Y = rpois(30, 3), data = X, family = "poisson"),
               "Family must be either")
})
