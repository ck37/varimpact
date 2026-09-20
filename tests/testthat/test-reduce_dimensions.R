# Tests for reduce_dimensions(), which was previously a one-line TODO.
#
# These lock in current behavior. Two of them are deliberately narrow, and the
# comments say why - see the notes on max_variables and newX below.
library(varimpact)
library(testthat)

context("reduce_dimensions")

make_df = function(p, n = 120L, seed = 11L) {
  set.seed(seed, "L'Ecuyer-CMRG")
  as.data.frame(matrix(rnorm(n * p), n, p))
}

test_that("data is returned untouched when it is already small enough", {
  data = make_df(5L)
  newX = make_df(5L, seed = 12L)

  result = reduce_dimensions(data, newX, max_variables = 10L)

  expect_identical(result$data, data)
  expect_equal(ncol(result$newX), ncol(newX))
  expect_identical(result$variables, colnames(data))
})

test_that("constant columns are dropped from both data and newX", {
  data = make_df(5L)
  newX = make_df(5L, seed = 12L)
  # A column with no variance carries no information and breaks downstream
  # models, so it should go - and it has to go from newX as well, otherwise
  # the two frames stop lining up.
  data$V3 = 1
  newX$V3 = 1

  result = reduce_dimensions(data, newX, max_variables = 10L)

  expect_false("V3" %in% colnames(result$data))
  expect_false("V3" %in% colnames(result$newX))
  expect_identical(colnames(result$data), c("V1", "V2", "V4", "V5"))
  expect_identical(result$variables, colnames(result$data))
})

test_that("dimensions are reduced when there are more columns than the cutoff", {
  data = make_df(20L)
  newX = make_df(20L, seed = 12L)

  result = reduce_dimensions(data, newX, max_variables = 5L)

  # Note: max_variables is not an upper bound. HOPACH picks the highest level
  # of the tree holding at least max_variables clusters, so the result can be
  # wider than max_variables - it is only guaranteed to be narrower than the
  # input. Asserting <= max_variables here would be asserting something the
  # function does not do.
  expect_lt(ncol(result$data), ncol(data))
  expect_gte(ncol(result$data), 1L)
  # Every returned column must be one of the originals, not a synthesized one.
  expect_true(all(colnames(result$data) %in% colnames(data)))
  expect_identical(result$variables, colnames(result$data))
})

test_that("newX stays aligned with data after reduction", {
  data = make_df(20L)
  newX = make_df(20L, seed = 12L)

  result = reduce_dimensions(data, newX, max_variables = 5L)

  # Callers feed both frames to the same fitted model, so the columns have to
  # match in name *and* order.
  expect_identical(colnames(result$data), colnames(result$newX))
  expect_equal(nrow(result$newX), nrow(newX))
})

test_that("columns missing from newX are added as zeros in data's order", {
  data = make_df(4L)
  # newX is missing V4 entirely, as happens when a factor level present in
  # training is absent from a validation fold.
  newX = make_df(4L, seed = 12L)[, 1:3]

  result = reduce_dimensions(data, newX, max_variables = 10L)

  expect_identical(colnames(result$newX), colnames(result$data))
  expect_true(all(result$newX$V4 == 0))
  # The columns that were already there keep their values.
  expect_equal(result$newX$V1, newX$V1)
})

test_that("a NULL max_variables disables reduction", {
  data = make_df(20L)
  newX = make_df(20L, seed = 12L)

  result = reduce_dimensions(data, newX, max_variables = NULL)

  expect_equal(ncol(result$data), ncol(data))
  expect_identical(result$variables, colnames(data))
})

test_that("verbose reports the reduction and the columns added to newX", {
  data = make_df(20L)
  newX = make_df(20L, seed = 12L)
  expect_output(reduce_dimensions(data, newX, max_variables = 5L, verbose = TRUE),
                "Reducing dimensions via clustering")

  data = make_df(4L)
  newX = make_df(4L, seed = 12L)[, 1:3]
  expect_output(reduce_dimensions(data, newX, max_variables = 10L, verbose = TRUE),
                "Adding missing columns in prediction data: V4")
})

test_that("when HOPACH fails on every attempt the data is returned unreduced", {
  # reduce_dimensions() retries hopach() with different settings, then gives
  # up and keeps the full data. hopach itself does not fail on well-formed
  # input, so make it.
  data = make_df(20L)
  newX = make_df(20L, seed = 12L)
  testthat::local_mocked_bindings(
    hopach = function(...) stop("synthetic HOPACH failure"),
    .package = "hopach")

  expect_output(
    result <- reduce_dimensions(data, newX, max_variables = 5L, verbose = TRUE),
    "Attempt 3 fail")
  expect_identical(result$data, data)
  expect_equal(ncol(result$newX), ncol(newX))
})

test_that("verbose reports constant columns being removed first", {
  data = make_df(5L)
  data$V2 = 3
  expect_output(reduce_dimensions(data, data, max_variables = 10L, verbose = TRUE),
                "First removing 1 constant columns")
})

test_that("newX really is optional", {
  # The signature has always said so, but the step that adds missing
  # columns to newX assumed there was one and failed on NULL.
  data = make_df(5L)
  data$V2 = 3
  result = reduce_dimensions(data, NULL, max_variables = 10L)
  expect_identical(colnames(result$data), c("V1", "V3", "V4", "V5"))
  expect_null(result$newX)
})

test_that("a failed distance matrix is reported and the data returned unreduced", {
  data = make_df(20L)
  newX = make_df(20L, seed = 12L)
  testthat::local_mocked_bindings(
    distancematrix = function(...) stop("synthetic distance failure"),
    .package = "hopach")
  expect_output(
    result <- reduce_dimensions(data, newX, max_variables = 5L, verbose = TRUE),
    "failed to calculate distance matrix")
  expect_identical(result$data, data)
})
