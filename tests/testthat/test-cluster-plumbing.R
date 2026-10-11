library(varimpact)
library(testthat)

context("cluster ids through validation, pooling and compile_results")

set.seed(31, "L'Ecuyer-CMRG")
n = 160
W = data.frame(w1 = rnorm(n), w2 = rnorm(n))
A = rbinom(n, 1, plogis(0.4 * W$w1))
Y = rbinom(n, 1, plogis(0.5 * W$w1 + A))
lib = c("SL.mean", "SL.glm")

fit = estimate_tmle2(Y = Y, A = A, W = W, family = "binomial",
                     Q.lib = lib, g.lib = lib, V = 2, Qbounds = c(0, 1))

test_that("apply_tmle_to_validation returns the row index it was given", {
  idx = sort(sample(n, 50))
  preds = varimpact:::apply_tmle_to_validation(Y[idx], A[idx], W[idx, ],
                                               "binomial", tmle = fit, id = idx)
  expect_identical(preds$row_index, idx)
  # The default is positional.
  preds = varimpact:::apply_tmle_to_validation(Y[idx], A[idx], W[idx, ],
                                               "binomial", tmle = fit)
  expect_identical(preds$row_index, seq_along(idx))
  expect_error(varimpact:::apply_tmle_to_validation(Y[idx], A[idx], W[idx, ],
                                                    "binomial", tmle = fit,
                                                    id = idx[-1]),
               "id has")
})

test_that("estimate_pooled_results returns rows aligned with each fold's curve", {
  # Shuffled halves, so row order within a fold is not 1..n.
  halves = split(sample(n), rep(1:2, length.out = n))
  fold_results = lapply(halves, function(idx) {
    list(val_preds = varimpact:::apply_tmle_to_validation(
      Y[idx], A[idx], W[idx, ], "binomial", tmle = fit, id = idx))
  })
  pooled = varimpact:::estimate_pooled_results(fold_results)
  expect_length(pooled$rows, 2L)
  for (k in 1:2) {
    expect_identical(pooled$rows[[k]], halves[[k]])
    expect_length(pooled$influence_curves[[k]], length(halves[[k]]))
  }

  # An empty fold keeps its slot empty, as for the curves.
  fold_results[[2]]$val_preds = NULL
  pooled = varimpact:::estimate_pooled_results(fold_results)
  expect_length(pooled$rows, 2L)
  expect_null(pooled$rows[[2]])
  expect_identical(pooled$rows[[1]], halves[[1]])
})

test_that("pooling without a row_index column still works", {
  idx = 1:80
  preds = varimpact:::apply_tmle_to_validation(Y[idx], A[idx], W[idx, ],
                                               "binomial", tmle = fit)
  preds$row_index = NULL
  pooled = varimpact:::estimate_pooled_results(list(list(val_preds = preds)))
  expect_null(pooled$rows)
  expect_length(pooled$influence_curves, 1L)
})

test_that("estimate_tmle2 accepts clustered ids, including very few clusters", {
  cl = rep(1:40, each = 4)
  expect_error(estimate_tmle2(Y = Y, A = A, W = W, family = "binomial",
                              Q.lib = lib, g.lib = lib, V = 10,
                              Qbounds = c(0, 1), id = cl), NA)
  # Fewer clusters than SuperLearner folds must not break CVFolds().
  few = rep(1:3, length.out = n)
  expect_error(estimate_tmle2(Y = Y, A = A, W = W, family = "binomial",
                              Q.lib = lib, g.lib = lib, V = 10,
                              Qbounds = c(0, 1), id = few), NA)
})

# Same fixture helper as test-compile-results.R.
make_vim = function(name, theta, var_ic, theta_rr, var_ic_log_rr, V = 2L) {
  list(EY1V = rep(0.5 + theta / 2, V),
       EY0V = rep(0.5 - theta / 2, V),
       thetaV = rep(theta, V),
       thetaV_rr = rep(theta_rr, V),
       varICV = rep(var_ic, V),
       varICV_log_rr = rep(var_ic_log_rr, V),
       labV = matrix(rep(c("lo", "hi"), each = V), nrow = V),
       nV = rep(50L, V),
       type = "factor",
       name = name)
}

vims = list(a = make_vim("a", 0.3, 1, 1.2, 1),
            b = make_vim("b", 0.1, 1, 1.05, 1))
compile = function(...) {
  varimpact:::compile_results(colnames_numeric = character(0),
                              colnames_factor = names(vims),
                              vim_numeric = list(),
                              vim_factor = unname(vims),
                              V = 2L, ...)
}

test_that("compile_results with df = Inf is unchanged", {
  expect_identical(compile(df = Inf), compile())
})

test_that("finite df uses t critical values and p-values", {
  df = 9
  raw = compile(df = df)$results_raw
  se = sqrt(1 / 100)
  expect_equal(raw$rawp, stats::pt(c(0.3, 0.1) / se, df, lower.tail = FALSE))
  expect_equal(raw$rr_rawp,
               stats::pt(log(c(1.2, 1.05)) / se, df, lower.tail = FALSE))
  # t intervals are wider than normal ones.
  norm_raw = compile()$results_raw
  expect_true(all(raw$rawp > norm_raw$rawp))
  expect_identical(raw$AvePsi, norm_raw$AvePsi)
  crit = stats::qt(0.975, df)
  expected_ci = paste0("(", signif(0.3 - crit * se, 3), " - ",
                       signif(0.3 + crit * se, 3), ")")
  expect_identical(raw["a", "CI95"], expected_ci)
})
