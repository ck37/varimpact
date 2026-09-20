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
