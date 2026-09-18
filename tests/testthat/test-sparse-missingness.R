library(varimpact)
library(testthat)

# A missingness indicator with only a handful of zeros cannot support a
# covariate-adjusted fit. estimate_tmle2() should fall back to the marginal
# rate rather than handing the degenerate outcome to SuperLearner, where
# learners error one by one and the surviving ensemble is the intercept anyway.

context("Sparse missingness mechanism")

test_that("marginal_fit predicts its probability for any number of rows", {
  fit = varimpact:::marginal_fit(0.25)
  expect_s3_class(fit, "varimpact_marginal")
  expect_equal(predict(fit, data.frame(a = 1:4)), rep(0.25, 4))
  # Called the way predict_g_delta() calls a non-SuperLearner fit.
  expect_equal(stats::predict(fit, newdata = data.frame(a = 1:3),
                              type = "response"),
               rep(0.25, 3))
  expect_error(varimpact:::marginal_fit(c(0.1, 0.2)))
  expect_error(varimpact:::marginal_fit(1.5))
})

test_that("tmle_estimate_g falls back to the marginal rate below min_cell_size", {
  set.seed(1)
  n = 60
  W = data.frame(W1 = rnorm(n), W2 = rnorm(n))
  A = rbinom(n, 1, 0.5)
  # Two observations missing out of 60.
  delta = c(0, 0, rep(1, n - 2))
  d = data.frame(delta, Z = 1, A, W)

  g = varimpact:::tmle_estimate_g(d = d, SL.library = c("SL.mean", "SL.glmnet"),
                                  V = 2, stratify = TRUE, outcome = "D",
                                  min_cell_size = 4)

  expect_equal(g$type, "marginal")
  expect_s3_class(g$model, "varimpact_marginal")
  expect_equal(colnames(g$g1W), c("Z0A0", "Z0A1", "Z1A0", "Z1A1"))
  expect_true(all(g$g1W == mean(delta)))

  # And the fit carries over to a validation fold at the same value.
  expect_equal(varimpact:::predict_g_delta(g$model, W), rep(mean(delta), n))
})

test_that("the fallback reports itself when verbose", {
  set.seed(1)
  n = 60
  W = data.frame(W1 = rnorm(n), W2 = rnorm(n))
  A = rbinom(n, 1, 0.5)
  delta = c(0, 0, rep(1, n - 2))
  d = data.frame(delta, Z = 1, A, W)

  output = capture.output(
    varimpact:::tmle_estimate_g(d = d, SL.library = "SL.mean", V = 2,
                                stratify = TRUE, outcome = "D",
                                min_cell_size = 4, verbose = TRUE))

  expect_true(any(grepl("marginal proportion", output)))
  # The message names the count that tripped the guard and the floor it missed.
  expect_true(any(grepl("2 observations, fewer than 4", output)))
})

test_that("min_cell_size = 0 leaves the covariate-adjusted fit alone", {
  set.seed(2)
  n = 60
  W = data.frame(W1 = rnorm(n), W2 = rnorm(n))
  A = rbinom(n, 1, 0.5)
  delta = rbinom(n, 1, 0.5)
  d = data.frame(delta, Z = 1, A, W)

  g = suppressWarnings(
    varimpact:::tmle_estimate_g(d = d, SL.library = "SL.mean", V = 2,
                                stratify = TRUE, outcome = "D"))

  expect_false(identical(g$type, "marginal"))
})

test_that("a sparsely missing covariate does not error out of glmnet", {
  skip_on_cran()
  set.seed(1, "L'Ecuyer-CMRG")
  future::plan("sequential")

  n = 200
  X = as.data.frame(matrix(rnorm(n * 4), n, 4))
  Y = rbinom(n, 1, plogis(0.2 * X[, 1] - 0.2 * X[, 3]))
  # A handful of missing covariate values, as in the README example. Each one
  # makes delta = 0 for that row once the variable is binned.
  for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

  # SuperLearner prints a learner's error rather than signalling it, so the
  # errors only show up on the message stream.
  output = capture.output(
    vim <- suppressWarnings(
      varimpact(Y = Y, data = X, Q.library = c("SL.mean", "SL.glmnet"),
                g.library = c("SL.mean", "SL.glmnet"), verbose = FALSE)),
    type = "message")

  expect_false(any(grepl("lognet", output)))
  expect_false(any(grepl("non-conformable", output)))
  expect_true(nrow(vim$results_all) > 0)
})

test_that("tmle_estimate_g() falls back to glm when SuperLearner fails", {
  # SuperLearner is wrapped in try(), so a library that cannot be resolved at
  # all leaves m a try-error rather than stopping the run. The function then
  # fits a main terms glm instead and carries on, which is what keeps one bad
  # library from ending an entire varimpact() call.
  #
  # This also exercises the non-SuperLearner arm of the counterfactual
  # reshaping: predict() is called without the X and Y arguments a
  # SuperLearner fit needs, and the four Z/A combinations still come back.
  set.seed(1)
  n <- 60
  W <- data.frame(W1 = rnorm(n), W2 = rnorm(n))
  A <- rbinom(n, 1, 0.5)
  delta <- rbinom(n, 1, 0.8)
  d <- data.frame(delta, Z = 1, A, W)

  output <- capture.output(
    g <- suppressWarnings(
      varimpact:::tmle_estimate_g(d = d,
                                  SL.library = "SL.thisLearnerDoesNotExist",
                                  V = 2, stratify = TRUE, outcome = "D",
                                  verbose = TRUE,
                                  message = "missingness mechanism")))

  # The failure announces itself unconditionally, not only when verbose; the
  # description of what replaced it is verbose-only.
  expect_true(any(grepl("Defaulting to glm", output, fixed = TRUE)))
  expect_true(any(grepl("Running main terms regression", output, fixed = TRUE)))

  expect_equal(g$type, "glm")
  expect_s3_class(g$model, "glm")
  expect_equal(colnames(g$g1W), c("Z0A0", "Z0A1", "Z1A0", "Z1A1"))
  expect_equal(dim(g$g1W), c(n, 4L))
  expect_true(all(g$g1W >= 0 & g$g1W <= 1))

  # coef comes from the glm rather than staying the scalar NA a SuperLearner
  # fit leaves it as, which is how a caller can tell which path ran.
  expect_named(g$coef, c("(Intercept)", "Z", "A", "W1", "W2"))
  expect_true(is.finite(g$coef[["(Intercept)"]]))

  # Z is NA, and that is the design rather than a defect: the missingness
  # mechanism is fit on data.frame(delta, Z = 1, A, W), so Z is constant and
  # the fit is rank deficient in it. This is the same rank deficiency behind
  # the "prediction from rank-deficient fit" warnings the suite emits, and
  # predict_g_delta() suppresses for the same reason.
  expect_true(is.na(g$coef[["Z"]]))
})
