library(varimpact)
library(testthat)

context("Imputation helpers")

# impute_median() and impute_knn() replaced caret::preProcess() with
# method = "medianImpute" and "knnImpute". Both were checked against caret
# 7.0-1 when written and matched it to the bit; the knn snapshot below records
# caret's values so that the equivalence stays under test without caret.

make_data = function(n = 120) {
  set.seed(3, "L'Ecuyer-CMRG")
  d = data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  d$a[c(2, 7, 30)] = NA
  d$b[c(5, 7)] = NA
  d$c[c(11, 40, 55)] = NA
  d
}

# The knn algorithm written out the slow way: center and scale, then for
# every incomplete row rank the complete rows by Euclidean distance over the
# columns the row does have, and average the k nearest in each column it lacks.
knn_reference = function(d, k = 5L) {
  x = as.matrix(d)
  x = sweep(x, 2, apply(x, 2, mean, na.rm = TRUE), "-")
  x = sweep(x, 2, apply(x, 2, sd, na.rm = TRUE), "/")
  complete = x[complete.cases(x), , drop = FALSE]
  for (i in which(!complete.cases(x))) {
    observed = !is.na(x[i, ])
    distances = sqrt(colSums((t(complete[, observed, drop = FALSE]) - x[i, observed])^2))
    nearest = complete[order(distances)[seq_len(k)], , drop = FALSE]
    x[i, !observed] = apply(nearest[, !observed, drop = FALSE], 2, mean)
  }
  x
}

test_that("median imputation fills each column with its own median", {
  d = make_data()
  out = varimpact:::impute_median(d)

  expect_s3_class(out$data, "data.frame")
  expect_identical(dim(out$data), dim(d))
  expect_identical(names(out$data), names(d))
  expect_false(anyNA(out$data))
  expect_equal(out$medians, vapply(d, median, numeric(1), na.rm = TRUE))
  # Observed values are untouched, and only they.
  observed = !is.na(d)
  expect_identical(as.matrix(out$data)[observed], as.matrix(d)[observed])
  for (col in names(d)) {
    expect_true(all(out$data[[col]][is.na(d[[col]])] == median(d[[col]], na.rm = TRUE)))
  }
})

test_that("median imputation keeps integer columns whole", {
  d = data.frame(i = c(1L, NA, 3L, 4L, 9L), x = c(NA, 0.5, 1.5, 2, 3))
  out = varimpact:::impute_median(d)
  expect_equal(out$data$i, c(1, 3.5, 3, 4, 9))
  expect_equal(out$data$x, c(1.75, 0.5, 1.5, 2, 3))
})

test_that("knn imputation matches a brute-force neighbor search", {
  d = make_data()
  out = varimpact:::impute_knn(d)

  expect_s3_class(out$data, "data.frame")
  expect_identical(dim(out$data), dim(d))
  expect_identical(names(out$data), names(d))
  expect_false(anyNA(out$data))
  expect_equal(out$k, 5L)
  expect_equal(as.matrix(out$data), knn_reference(d), check.attributes = FALSE)
})

test_that("knn output is centered and scaled", {
  d = make_data()
  out = varimpact:::impute_knn(d)
  expect_equal(out$center, vapply(d, mean, numeric(1), na.rm = TRUE))
  expect_equal(out$scale, vapply(d, sd, numeric(1), na.rm = TRUE))
  scaled = as.matrix(out$data)
  observed = !is.na(d)
  expected = sweep(sweep(as.matrix(d), 2, out$center, "-"), 2, out$scale, "/")
  expect_equal(scaled[observed], expected[observed])
})

test_that("knn imputation reproduces caret's values", {
  # Recorded from predict(caret::preProcess(d, method = "knnImpute"), d) with
  # caret 7.0-1 and RANN 2.6.2, for the eight missing cells of make_data().
  d = make_data()
  out = as.matrix(varimpact:::impute_knn(d)$data)
  cells = which(is.na(d), arr.ind = TRUE)
  caret_values = c(0.08242208095635256, -0.00070744258585666794, 0.28595067529422985,
                   0.33875807714007494, 1.4670421574577595,
                   -0.50017613452599907, 1.2174898642583405, 0.53505829072519839)
  expect_equal(out[cells], caret_values, tolerance = 1e-12)
})

test_that("knn treats a constant column as centered only", {
  d = make_data()
  d$flat = 1
  # Row 7 is already incomplete, so the set of complete rows that neighbors
  # are drawn from is the same with or without the flat column.
  d$flat[7] = NA
  out = varimpact:::impute_knn(d)
  expect_equal(unname(out$scale["flat"]), 1)
  expect_true(all(out$data$flat == 0))
  # The other columns come out as they would without the flat column: its
  # distances are all zero, so it changes no neighbor.
  expect_equal(as.matrix(out$data[, c("a", "b", "c")]),
               knn_reference(d[, c("a", "b", "c")]), check.attributes = FALSE)
})

test_that("knn uses fewer neighbors when fewer complete rows exist", {
  # caret errored here (RANN cannot return more neighbors than points).
  d = data.frame(a = c(1, 2, 3, NA, 5, 6), b = c(2, 4, 6, 8, NA, 12))
  out = varimpact:::impute_knn(d)
  expect_equal(out$k, 4L)
  expect_false(anyNA(out$data))
  expect_equal(as.matrix(out$data), knn_reference(d, k = 4L), check.attributes = FALSE)

  none_complete = data.frame(a = c(1, NA, 3), b = c(NA, 2, NA))
  expect_error(varimpact:::impute_knn(none_complete), "no row is complete")
})
