# Test single column support in varimpact
# This addresses issue #6: "Single column in data"

library(testthat)
library(varimpact)

context("Single column support")

# Helper: run an expression, collecting warning messages rather than emitting
# them. varimpact() prints progress even with verbose = FALSE, and the
# estimation can warn about rank-deficient fits, so expect_silent() can never
# hold here; we assert on the specific warning we care about instead.
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cond) {
      warnings <<- c(warnings, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = warnings)
}

single_var_warning <- "Using single variable for variable importance analysis"

test_that("separate_factors_numerics handles single columns", {
  # Test single numeric column
  data_numeric <- data.frame(x1 = c(1.1, 2.2, 3.3, 4.4, 5.5))
  result <- separate_factors_numerics(data_numeric)
  expect_equal(ncol(result$df_factors), 0)
  expect_equal(ncol(result$df_numerics), 1)
  expect_equal(colnames(result$df_numerics), "x1")

  # Test single factor column
  data_factor <- data.frame(x1 = factor(c("A", "B", "C", "A", "B")))
  result <- separate_factors_numerics(data_factor)
  expect_equal(ncol(result$df_factors), 1)
  expect_equal(ncol(result$df_numerics), 0)
  expect_equal(colnames(result$df_factors), "x1")

  # Test single character column (should be converted to factor)
  data_char <- data.frame(x1 = c("A", "B", "C", "A", "B"), stringsAsFactors = FALSE)
  result <- separate_factors_numerics(data_char)
  expect_equal(ncol(result$df_factors), 1)
  expect_equal(ncol(result$df_numerics), 0)
  expect_equal(colnames(result$df_factors), "x1")
  expect_true(is.factor(result$df_factors$x1))
})

test_that("varimpact handles vector input", {
  # Create test data. N needs to be large enough that both CV folds retain
  # more than one bin after penalized histogramming, otherwise the variable is
  # legitimately skipped and results_all is NULL.
  set.seed(1)
  N <- 200
  X_vector <- rnorm(N)
  Y <- rbinom(N, 1, plogis(0.2 * X_vector))

  # Test that vector input is converted to data frame and processed
  expect_warning(
    vim <- varimpact(Y = Y, data = X_vector, verbose = FALSE, V = 2L),
    single_var_warning
  )

  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_s3_class(vim$results_all, "data.frame")
  expect_equal(nrow(vim$results_all), 1)
})

test_that("verbose reports the conversion and the single variable", {
  # varimpact()'s argument checks have two verbose-only branches: one
  # announcing that a vector was converted to a data frame, the other that only
  # one variable is present. Both are worth saying out loud, because each
  # changes what the run is actually doing relative to what was passed in.
  set.seed(1)
  N <- 200
  X_vector <- rnorm(N)
  Y <- rbinom(N, 1, plogis(0.2 * X_vector))

  output <- capture.output(
    expect_warning(
      varimpact(Y = Y, data = X_vector, verbose = TRUE, V = 2L),
      single_var_warning
    ))

  expect_true(any(grepl("Converting vector input to data frame", output,
                        fixed = TRUE)))
  expect_true(any(grepl("Single variable detected in data", output,
                        fixed = TRUE)))
})

test_that("varimpact handles single column data frame", {
  # Create test data
  set.seed(1)
  N <- 200
  X_single <- data.frame(x1 = rnorm(N))
  Y <- rbinom(N, 1, plogis(0.2 * X_single$x1))

  # Test that single column data frame is processed with warning
  expect_warning(
    vim <- varimpact(Y = Y, data = X_single, verbose = FALSE, V = 2L),
    single_var_warning
  )

  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_s3_class(vim$results_all, "data.frame")
  expect_equal(nrow(vim$results_all), 1)
  expect_equal(rownames(vim$results_all), "x1")
})

test_that("varimpact handles single factor column", {
  # Create test data
  set.seed(1)
  N <- 200
  X_factor <- data.frame(x1 = factor(sample(c("A", "B", "C"), N, replace = TRUE)))
  Y <- rbinom(N, 1, plogis(0.2 * as.numeric(X_factor$x1)))

  # Test that single factor column is processed with warning
  expect_warning(
    vim <- varimpact(Y = Y, data = X_factor, verbose = FALSE, V = 2L),
    single_var_warning
  )

  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_s3_class(vim$results_all, "data.frame")
  expect_equal(nrow(vim$results_all), 1)
})

test_that("varimpact still works with multiple columns (regression test)", {
  # Create test data
  set.seed(1)
  N <- 200
  X_multi <- data.frame(
    x1 = rnorm(N),
    x2 = rnorm(N),
    x3 = factor(sample(c("A", "B"), N, replace = TRUE))
  )
  Y <- rbinom(N, 1, plogis(0.2 * X_multi$x1 + 0.1 * X_multi$x2))

  # Multiple columns must not trigger the single-variable warning.
  run <- collect_warnings(
    vim <- varimpact(Y = Y, data = X_multi, verbose = FALSE, V = 2L)
  )
  expect_false(any(grepl(single_var_warning, run$warnings, fixed = TRUE)))

  # Check that the result is valid
  expect_s3_class(vim, "varimpact")
  expect_s3_class(vim$results_all, "data.frame")
  expect_equal(nrow(vim$results_all), 3)
})

test_that("tmle_estimate_g() falls back to SL.mean with no adjustment variables", {
  # When the analyzed variable is the only column in the data there is nothing
  # left to adjust for, so d reaches tmle_estimate_g() with a single column and
  # the covariate matrix is zero-column. Learners that regress on the covariates
  # cannot fit that at all - SL.glm's "Y ~ ." errors with "'.' in formula and no
  # 'data' argument" - so the library is narrowed to SL.mean up front rather
  # than letting each learner fail and be dropped.
  set.seed(1)
  d <- data.frame(A = rbinom(40, 1, 0.5))

  output <- capture.output(
    g <- varimpact:::tmle_estimate_g(d = d,
                                     SL.library = c("SL.mean", "SL.glm"),
                                     V = 2, stratify = TRUE, verbose = TRUE,
                                     message = "treatment mechanism"))

  expect_true(any(grepl("No adjustment variables", output, fixed = TRUE)))
  expect_true(any(grepl("with SL.mean alone", output, fixed = TRUE)))

  # SL.glm is in the library and would have errored; the fit succeeds anyway.
  expect_equal(g$type, "SuperLearner")

  # The conditional probability reduces to the marginal one, so every fitted
  # value is the sample mean. That is the same answer the dropped learners
  # would have left behind, arrived at without the failures.
  expect_equal(unique(as.vector(g$g1W)), mean(d$A))

  # Quiet by default.
  expect_silent(
    varimpact:::tmle_estimate_g(d = d, SL.library = "SL.mean", V = 2,
                                stratify = TRUE, message = "treatment mechanism"))
})
