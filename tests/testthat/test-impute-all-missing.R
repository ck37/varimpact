library(varimpact)
library(testthat)

context("Imputation when a row is missing every numeric covariate")

# A row in which every column is missing has no observed column to match
# neighbors on. The caret-based knn imputation this package used to call
# stopped on such a row, which took down the whole run for impute = "knn",
# while median and zero imputation filled the same row and carried on
# (issue #7).

make_data = function(n = 120, all_missing_rows = integer(0)) {
  set.seed(3, "L'Ecuyer-CMRG")
  d = data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  # Ordinary partial missingness, present in every case.
  d$a[c(2, 7, 30)] = NA
  d$b[c(5, 7)] = NA
  d$c[c(11, 40, 55)] = NA
  for (i in all_missing_rows) d[i, ] = NA
  d
}

impute_with = function(d, method, verbose = FALSE) {
  varimpact:::process_numerics(d,
                               quantile_probs_numeric = NULL,
                               miss.cut = 0.5,
                               bins_numeric = 10L,
                               impute = method,
                               verbose = verbose)$data.numW
}

test_that("every imputation method accepts an all-missing row", {
  d = make_data(all_missing_rows = c(13, 77))
  for (method in c("knn", "median", "zero")) {
    out = impute_with(d, method)
    expect_false(anyNA(out), info = method)
    expect_equal(nrow(out), nrow(d), info = method)
  }
})

test_that("knn puts all-missing rows at the column mean", {
  d = make_data(all_missing_rows = c(13, 77))
  all_missing = rowSums(is.na(d)) == ncol(d)
  out = as.matrix(impute_with(d, "knn"))
  # knn imputation centers and scales, so the column mean is 0.
  expect_true(all(out[all_missing, ] == 0))
})

test_that("all-missing rows do not disturb the other rows", {
  # The rows that can be imputed must come out exactly as they would from a
  # frame that never had the all-missing rows - those are not complete cases,
  # so they are never anyone's neighbor, and being NA in every column they
  # contribute to no column's mean or sd either.
  d = make_data(all_missing_rows = c(13, 77))
  all_missing = rowSums(is.na(d)) == ncol(d)
  out = as.matrix(impute_with(d, "knn"))

  sub = d[!all_missing, , drop = FALSE]
  reference = as.matrix(impute_with(sub, "knn"))
  expect_equal(out[!all_missing, ], reference,
               tolerance = 0, check.attributes = FALSE)
})

test_that("process_numerics() returns what impute_knn() computes", {
  # Guards the common case: with nothing to skip, the numeric processing must
  # hand back the imputation helper's output as is, to the bit.
  d = make_data()
  expect_false(any(rowSums(is.na(d)) == ncol(d)))
  out = as.matrix(impute_with(d, "knn"))
  reference = as.matrix(varimpact:::impute_knn(d)$data)
  expect_equal(out, reference, tolerance = 0, check.attributes = FALSE)
})

test_that("knn reports the rows it could not impute from neighbors", {
  # These rows are filled with a value no neighbor supplied, so saying so when
  # the caller asked for output is part of the behaviour, not decoration.
  d = make_data(all_missing_rows = c(13, 77))
  expect_output(impute_with(d, "knn", verbose = TRUE),
                "2 row\\(s\\) are missing every numeric covariate")

  # And it stays quiet when there is nothing to report.
  expect_silent(impute_with(make_data(), "knn", verbose = FALSE))
})

test_that("varimpact() itself runs with an all-missing row", {
  # The end-to-end path issue #7 reported against.
  set.seed(1, "L'Ecuyer-CMRG")
  future::plan("sequential")
  n = 200
  X = data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  Y = rbinom(n, 1, plogis(X$x1))
  X[5, ] = NA                     # every covariate missing
  X$x1[9] = NA; X$x2[9] = NA      # all but one missing

  vim = varimpact(Y = Y, data = X, V = 2L, impute = "knn",
                  Q.library = c("SL.mean", "SL.glm"),
                  g.library = c("SL.mean", "SL.glm"),
                  verbose = FALSE)
  expect_s3_class(vim, "varimpact")
  expect_true(nrow(vim$results_all) > 0)
})
