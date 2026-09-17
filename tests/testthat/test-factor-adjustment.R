# Regression tests for the adjustment matrix used when estimating the importance
# of a numeric variable.
#
# vim_numerics() builds W from factors$datafac.dumW, but process_factors() only
# returned datafac.dum - there was no datafac.dumW key. The expression was
# therefore always NULL, was replaced by rep(NA, n.cont), and was then dropped by
# the all-NA filter, so factor covariates were silently never adjusted for when
# estimating a numeric variable's importance.

library(testthat)
library(varimpact)

context("Factor covariates in the adjustment set")

make_factors <- function(df) {
  varimpact:::process_factors(varimpact:::separate_factors_numerics(df)$df_factors,
                              quantile_probs_factor = c(0.1, 0.9),
                              miss.cut = 0.5,
                              verbose = FALSE)
}

test_that("process_factors() returns the imputed indicator matrix", {
  set.seed(1)
  df <- data.frame(F1 = factor(sample(c("a", "b", "c"), 50, replace = TRUE)))
  f <- make_factors(df)

  # vim_numerics() reads this key; without it the factor columns never reach W.
  expect_true("datafac.dumW" %in% names(f))
  expect_false(is.null(f$datafac.dumW))
  expect_equal(colnames(f$datafac.dumW), colnames(f$datafac.dum))
  expect_equal(nrow(f$datafac.dumW), nrow(df))
})

test_that("the adjustment matrix carries no missing values", {
  set.seed(2)
  df <- data.frame(F1 = factor(sample(c("a", "b", "c"), 50, replace = TRUE)))
  df$F1[c(3, 17, 40)] <- NA
  f <- make_factors(df)

  # factors_to_indicators() encodes a missing value as its own level rather than
  # leaving NA in the dummies, so datafac.dum is already complete here. The
  # imputation in datafac.dumW is a guarantee for the paths where it is not;
  # what matters downstream is that nothing NA reaches the adjustment set.
  expect_true("F1XXNA" %in% colnames(f$datafac.dum))
  expect_true("Imiss_F1" %in% colnames(f$miss.fac))
  expect_false(any(is.na(f$datafac.dumW)))
  expect_equal(f$datafac.dumW[!is.na(f$datafac.dum)],
               f$datafac.dum[!is.na(f$datafac.dum)])
})

test_that("process_factors() returns NULL for both when there are no factors", {
  f <- make_factors(data.frame(V1 = rnorm(20)))
  expect_null(f$datafac.dum)
  expect_null(f$datafac.dumW)
  expect_equal(f$num_factors, 0L)
})

test_that("a numeric variable is adjusted for factor covariates", {
  # The end-to-end consequence. Traced through reduce_dimensions(), which
  # receives the adjustment matrix for each numeric candidate variable; before
  # the fix it saw only the other numerics.
  skip_if_not_installed("future")
  future::plan("sequential")

  # The tracer is evaluated inside the traced function, and varimpact() farms the
  # per-variable work out through future_lapply, so the accumulator has to live
  # somewhere both can reach: refer to the global environment explicitly rather
  # than relying on lexical scoping from this test frame.
  assign("vi_adjustment_cols", character(), envir = globalenv())
  on.exit(suppressWarnings(rm("vi_adjustment_cols", envir = globalenv())),
          add = TRUE)
  trace(varimpact:::reduce_dimensions, print = FALSE,
        tracer = quote({
          assign("vi_adjustment_cols",
                 unique(c(get("vi_adjustment_cols", envir = globalenv()),
                          colnames(data))),
                 envir = globalenv())
        }))
  on.exit(try(untrace(varimpact:::reduce_dimensions), silent = TRUE), add = TRUE)

  set.seed(1, "L'Ecuyer-CMRG")
  N <- 200
  X <- data.frame(V1 = rnorm(N), V2 = rnorm(N),
                  F1 = factor(sample(c("a", "b", "c"), N, replace = TRUE)))
  Y <- rbinom(N, 1, plogis(0.3 * X$V1))

  invisible(suppressWarnings(
    varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"))))

  untrace(varimpact:::reduce_dimensions)
  cols <- get("vi_adjustment_cols", envir = globalenv())

  expect_true(any(grepl("^F1XX", cols)),
              info = paste("adjustment columns seen:",
                           paste(sort(cols), collapse = ", ")))
  expect_true(all(c("V1", "V2") %in% cols))
})
