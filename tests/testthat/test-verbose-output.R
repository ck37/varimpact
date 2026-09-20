library(varimpact)
library(testthat)

# Every stage of varimpact() has verbose = TRUE reporting that no other test
# turns on. One small run with a numeric and a factor variable, and one with a
# factor that gets dropped, exercise that reporting end to end. The checks are
# for the messages that mark each stage, so a cat() with a broken argument
# would surface here rather than in a user's console.

context("verbose output")

future::plan("sequential")

set.seed(13, "L'Ecuyer-CMRG")
n = 120
X = data.frame(x1 = rnorm(n),
               x2 = rnorm(n),
               f1 = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
X$x2[1:4] = NA
Y = rbinom(n, 1, plogis(0.8 * X$x1))

test_that("a verbose run reports each stage", {
  out = capture.output(
    vim <- varimpact(Y = Y, data = X, V = 2L, verbose = TRUE,
                     Q.library = c("SL.mean", "SL.glm"),
                     g.library = c("SL.mean", "SL.glm"),
                     bins_numeric = 3L))
  expect_s3_class(vim, "varimpact")
  text = paste(out, collapse = "\n")
  for (marker in c("Processing numerics", "Processing factors", "Min level prediction",
                   "Completed numeric variable importance estimation")) {
    expect_true(grepl(marker, text, fixed = TRUE), info = marker)
  }
})

test_that("a verbose run says when every factor is dropped", {
  X_flat = data.frame(x1 = rnorm(n), x2 = rnorm(n), f = factor(rep("a", n)))
  out = capture.output(
    vim <- varimpact(Y = Y, data = X_flat, V = 2L, verbose = TRUE,
                     Q.library = "SL.mean", g.library = "SL.mean",
                     bins_numeric = 3L))
  expect_s3_class(vim, "varimpact")
  expect_true(any(grepl("All factors were dropped", out, fixed = TRUE)))
})

test_that("a verbose run reports variables left out of A_names and correlated columns", {
  set.seed(14, "L'Ecuyer-CMRG")
  n = 120
  f1 = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  X = data.frame(x1 = rnorm(n),
                 f2 = factor(sample(c("p", "q"), n, replace = TRUE)),
                 # Nearly an indicator of f1 == "c", so it is dropped from
                 # f1's adjustment set under the default corthres = 0.8. (The
                 # check correlates against f1's dummy columns, which omit
                 # the reference level "a", so the level has to be another.)
                 x_corr = as.numeric(f1 == "c") + rnorm(n, sd = 0.05),
                 f1 = f1)
  Y = rbinom(n, 1, plogis(0.5 * X$x1))
  out = capture.output(
    vim <- varimpact(Y = Y, data = X, V = 2L, verbose = TRUE,
                     A_names = c("x1", "f1"),
                     Q.library = "SL.mean", g.library = "SL.mean",
                     bins_numeric = 3L))
  text = paste(out, collapse = "\n")
  expect_s3_class(vim, "varimpact")
  expect_true(grepl("Skipping x_corr  as it is not in A_names", text, fixed = TRUE))
  expect_true(grepl("Skipping f2 as it is not in A_names", text, fixed = TRUE))
  expect_true(grepl("columns based on correlation threshold", text, fixed = TRUE))
})
