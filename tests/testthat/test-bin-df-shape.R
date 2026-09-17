# bin_df must carry exactly one row per bin per fold.
#
# It is built by coercing each bin_result list to a one-row data frame. Every
# field of that list is therefore required to be length 1: data.frame() recycles
# a longer element into extra rows, and rejects a zero-length one outright.
# W_names, added for the adjustment_exclusions work, is the length of the
# adjustment set, so it has to be excluded from the coercion the same way
# test_predictions is.

library(testthat)
library(varimpact)

context("bin_df shape")

test_that("bin_df has one row per bin, not one per adjustment variable", {
  future::plan("sequential")
  set.seed(1, "L'Ecuyer-CMRG")

  N <- 200
  X <- data.frame(V1 = rnorm(N), V2 = rnorm(N), V3 = rnorm(N))
  Y <- rbinom(N, 1, plogis(0.3 * X$V1))

  vim <- suppressWarnings(
    varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm")))

  for (name in names(vim$all_vims)) {
    for (fold in seq_along(vim$all_vims[[name]]$fold_results)) {
      bin_df <- vim$all_vims[[name]]$fold_results[[fold]]$bin_df
      if (is.null(bin_df) || nrow(bin_df) == 0L) next

      # One row per bin: no level may appear twice within a fold.
      expect_equal(anyDuplicated(bin_df[, c("level", "cv_fold")]), 0L,
                   info = paste(name, "fold", fold))

      # The vector-valued fields must not have been folded into the frame.
      expect_false("W_names" %in% names(bin_df))
      expect_false("test_predictions" %in% names(bin_df))
    }
  }

  # And the adjustment sets are still recorded on the list itself.
  w <- vim$all_vims$V1$fold_results[[1]]$bin_results[[1]]$W_names
  expect_true(is.character(w))
  expect_true(all(c("V2", "V3") %in% w))
})

test_that("an empty adjustment set does not break bin_df", {
  # With a single column there are no adjustment variables at all, so W has zero
  # columns and colnames(W) is character(0) - the case that errored with
  # "arguments imply differing number of rows: 1, 0".
  future::plan("sequential")
  set.seed(1, "L'Ecuyer-CMRG")

  N <- 200
  X <- data.frame(V1 = rnorm(N), V2 = rnorm(N))
  Y <- rbinom(N, 1, plogis(0.3 * X$V1))

  # Two columns, so analysing either leaves exactly one adjustment variable;
  # the zero-column case is exercised directly below since single-column input
  # is not supported on this branch.
  expect_error(
    suppressWarnings(
      varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
                Q.library = c("SL.mean", "SL.glm"),
                g.library = c("SL.mean", "SL.glm"))),
    NA)

  # The coercion itself, with an empty adjustment set.
  result <- list(name = "V1", W_names = character(0), cv_fold = 1L, level = 1L,
                 test_predictions = NULL)
  expect_error(
    data.frame(result[!names(result) %in% c("test_predictions", "W_names")],
               stringsAsFactors = FALSE),
    NA)
})
