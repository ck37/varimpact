library(varimpact)
library(testthat)

# results_by_level() averages the per-fold results down to one row per level.
# It must do that WITHOUT renaming any column: it selects cv_fold,
# train_cell_size and test_cell_size away immediately afterwards, and callers
# read test_theta_tmle and friends by name.
#
# This is the regression that commit aa20337 introduced and 1ee9675 reverted.
# summarize_all(list(mean = mean)) passes a *named* list, which makes dplyr
# append the name to every output column, so cv_fold became cv_fold_mean and
# the following select() could not find it. Nothing asserted the column names,
# so it reached master. These tests assert them.

context("results_by_level")

make_folds = function() {
  # Two CV folds for one variable with two levels, shaped like the frame
  # vim_numerics() builds.
  data.frame(
    name = rep("V1", 4),
    level = rep(c(1, 2), each = 2),
    level_label = rep(c("[0,1]", "(1,2]"), each = 2),
    cv_fold = rep(c(1, 2), 2),
    train_cell_size = c(10, 12, 14, 16),
    test_cell_size = c(5, 6, 7, 8),
    train_theta_tmle = c(0.1, 0.3, 0.5, 0.7),
    test_theta_tmle = c(0.2, 0.4, 0.6, 0.8),
    test_var_tmle = c(0.01, 0.03, 0.05, 0.07),
    train_msg = rep("ok", 4),
    test_msg = rep("ok", 4),
    stringsAsFactors = FALSE
  )
}

test_that("column names survive the aggregation", {
  res = varimpact:::results_by_level(make_folds())

  expect_s3_class(res, "data.frame")
  # Exactly the input columns, minus the ones selected away. No suffixes.
  expect_equal(sort(names(res)),
               sort(c("name", "level", "level_label",
                      "train_theta_tmle", "test_theta_tmle", "test_var_tmle")))
  expect_false(any(grepl("_mean$", names(res))))
  # The columns the following select() drops are gone, not renamed.
  expect_false(any(c("cv_fold", "train_cell_size", "test_cell_size",
                     "cv_fold_mean", "train_cell_size_mean",
                     "test_cell_size_mean") %in% names(res)))
})

test_that("one row per level, averaged across folds", {
  res = varimpact:::results_by_level(make_folds())

  expect_equal(nrow(res), 2L)
  res = res[order(res$level), ]
  expect_equal(res$level, c(1, 2))
  expect_equal(res$level_label, c("[0,1]", "(1,2]"))
  expect_equal(res$train_theta_tmle, c(0.2, 0.6))
  expect_equal(res$test_theta_tmle, c(0.3, 0.7))
  expect_equal(res$test_var_tmle, c(0.02, 0.06))
})

test_that("aggregating does not warn", {
  # The deprecated dplyr::funs() emitted a deprecation warning on every call,
  # which is where varimpact's wall of identical backtraces came from.
  expect_no_warning(varimpact:::results_by_level(make_folds()))
})

test_that("a failure is reported rather than signalled", {
  # The function wraps its body in tryCatch and returns NULL on error; plot_var
  # relies on that. Missing columns is the easiest way to trip it.
  expect_null(suppressWarnings(
    varimpact:::results_by_level(data.frame(name = "V1", level = 1))))
})
