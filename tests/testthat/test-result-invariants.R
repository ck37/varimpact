# Structural invariants on the varimpact() result object.
#
# These deliberately assert relationships that must hold for any dataset,
# rather than specific estimates. varimpact() is data-adaptive - fold
# assignment is random and SuperLearner cross-validates internally - so
# pinning exact numbers produces tests that break on an unrelated R or
# SuperLearner upgrade, which is how a suite ends up disabled instead of
# fixed.
library(varimpact)
library(testthat)

context("Result object invariants")

future::plan("sequential")

# Parse the "(lower - upper)" CI strings that compile_results() formats.
parse_ci = function(x) {
  matches = regmatches(x, regexec("^\\(([-0-9.e+]+) - ([-0-9.e+]+)\\)$", x))
  t(vapply(matches, function(parts) as.numeric(parts[2:3]), numeric(2)))
}

fit_once = function(seed = 21L, N = 250L, p = 4L, V = 2L) {
  set.seed(seed, "L'Ecuyer-CMRG")
  X = as.data.frame(matrix(rnorm(N * p), N, p))
  Y = rbinom(N, 1, plogis(0.6 * X[, 1]))
  # varimpact() is chatty; the output is not what we are testing.
  utils::capture.output(vim <- varimpact(Y = Y, data = X, V = V))
  list(vim = vim, X = X, Y = Y, V = V)
}

test_that("every reported row corresponds to a real input variable", {
  fit = fit_once()
  results = fit$vim$results_all

  expect_s3_class(fit$vim, "varimpact")
  expect_true(nrow(results) > 0L)
  # Variables can legitimately be dropped (minCell, minYs, no variation in a
  # fold), so this is a subset relationship rather than an equality.
  expect_true(all(rownames(results) %in% colnames(fit$X)))
  expect_lte(nrow(results), ncol(fit$X))
  expect_false(any(duplicated(rownames(results))))
})

test_that("a reported result carries no missing estimates or p-values", {
  results = fit_once()$vim$results_all

  # A row that made it into results_all has passed compile_results()' own
  # filtering, so an NA here means something failed silently upstream.
  expect_false(any(is.na(results$Estimate)))
  expect_false(any(is.na(results[["P-value"]])))
  expect_false(any(is.na(results[["Adj. p-value"]])))
  expect_true(all(is.finite(results$Estimate)))
})

test_that("p-values are probabilities and adjustment never decreases them", {
  results = fit_once()$vim$results_all

  expect_true(all(results[["P-value"]] >= 0 & results[["P-value"]] <= 1))
  expect_true(all(results[["Adj. p-value"]] >= 0 & results[["Adj. p-value"]] <= 1))
  # Holm/BH adjustment can only move a p-value up.
  expect_true(all(results[["Adj. p-value"]] >= results[["P-value"]] - 1e-12))
})

test_that("confidence intervals bracket their own point estimate", {
  results = fit_once()$vim$results_all

  ci = parse_ci(results$CI95)
  expect_false(any(is.na(ci)))
  expect_true(all(ci[, 1] <= results$Estimate))
  expect_true(all(results$Estimate <= ci[, 2]))

  # The risk ratio is a ratio of means, so it and its interval are positive,
  # and the interval brackets the estimate on the log scale it was built on.
  expect_true(all(results[["Est. RR"]] > 0))
  ci_rr = parse_ci(results[["CI95 RR"]])
  expect_true(all(ci_rr > 0))
  expect_true(all(ci_rr[, 1] <= results[["Est. RR"]]))
  expect_true(all(results[["Est. RR"]] <= ci_rr[, 2]))
})

test_that("the reported estimate is the mean of its per-fold estimates", {
  fit = fit_once()
  results = fit$vim$results_all
  by_fold = fit$vim$results_by_fold

  expect_true(setequal(rownames(by_fold), rownames(results)))

  fold_cols = grep("^Est_v", names(by_fold), value = TRUE)
  expect_equal(length(fold_cols), fit$V)

  # compile_results() defines the reported estimate as the mean across folds,
  # so this has to hold exactly. It is the single assertion most likely to
  # catch a fold being dropped or double-counted.
  fold_means = rowMeans(by_fold[rownames(results), fold_cols, drop = FALSE])
  expect_equal(unname(fold_means), results$Estimate)
})

test_that("cross-validation folds cover every observation exactly once", {
  fit = fit_once()

  expect_equal(length(fit$vim$cv_folds), length(fit$Y))
  expect_false(any(is.na(fit$vim$cv_folds)))
  expect_setequal(unique(fit$vim$cv_folds), seq_len(fit$V))
  # Stratified assignment should keep the folds close to balanced.
  fold_sizes = table(fit$vim$cv_folds)
  expect_lte(max(fold_sizes) - min(fold_sizes), 2L)
})

test_that("results_consistent is a subset of results_all", {
  fit = fit_once()
  results = fit$vim$results_all
  consistent = fit$vim$results_consistent

  expect_true(all(rownames(consistent) %in% rownames(results)))
  expect_lte(nrow(consistent), nrow(results))
  expect_identical(names(consistent), setdiff(names(results), "Consistent"))
})

test_that("V is honored", {
  fit = fit_once(V = 3L)

  expect_equal(fit$vim$V, 3L)
  expect_equal(length(grep("^Est_v", names(fit$vim$results_by_fold))), 3L)
  expect_setequal(unique(fit$vim$cv_folds), 1:3)
})
