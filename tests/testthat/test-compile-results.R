library(varimpact)
library(testthat)

context("compile_results")

# One per-variable result list, shaped like what vim_factors() hands to
# compile_results(): two folds, constant across folds so that the expected
# p-values are easy to compute by hand.
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

# compile_results() defines the one-sided p-value from the mean estimate and
# the mean influence-curve variance over n = the first variable's total
# validation size, 100 here.
p_rd = function(theta, var_ic) 1 - pnorm(theta / sqrt(var_ic / 100))
p_rr = function(theta_rr, var_ic_log_rr) 1 - pnorm(log(theta_rr) / sqrt(var_ic_log_rr / 100))

test_that("relative-risk p-values stay with their own variable", {
  # The risk-difference ranking is a, b, c and the relative-risk ranking is
  # the reverse, c, b, a. The rows are ordered by the risk-difference p-value,
  # and each row's relative-risk p-values must be that variable's own. The
  # multtest-based version pasted the relative-risk p-values in sorted by their
  # own order, so on this input a's row would have carried c's values.
  vims = list(a = make_vim("a", theta = 0.3, var_ic = 1, theta_rr = 1.1, var_ic_log_rr = 1),
              b = make_vim("b", theta = 0.2, var_ic = 1, theta_rr = 1.3, var_ic_log_rr = 1),
              c = make_vim("c", theta = 0.1, var_ic = 1, theta_rr = 1.5, var_ic_log_rr = 1))

  out = varimpact:::compile_results(colnames_numeric = character(0),
                                    colnames_factor = names(vims),
                                    vim_numeric = list(),
                                    vim_factor = unname(vims),
                                    V = 2L)
  raw = out$results_raw
  expect_identical(rownames(raw), c("a", "b", "c"))

  expected_rd = c(a = p_rd(0.3, 1), b = p_rd(0.2, 1), c = p_rd(0.1, 1))
  expected_rr = c(a = p_rr(1.1, 1), b = p_rr(1.3, 1), c = p_rr(1.5, 1))
  expect_equal(raw$rawp, unname(expected_rd))
  expect_equal(raw$rr_rawp, unname(expected_rr))
  expect_equal(raw$AvePsi_rr, c(1.1, 1.3, 1.5))

  # The adjusted columns are the adjustment of each variable's own p-value.
  expect_equal(raw$rr_Holm, unname(p.adjust(expected_rr, "holm")))
  expect_equal(raw$rr_BH, unname(p.adjust(expected_rr, "BH")))
  expect_equal(raw$Holm, unname(p.adjust(expected_rd, "holm")))
  expect_equal(raw$BH, unname(p.adjust(expected_rd, "BH")))

  # And the same in the user-facing table.
  all = out$results_all
  expect_equal(all[["P-value RR"]], unname(expected_rr[rownames(all)]))
  expect_equal(all[["Adj. p-value RR"]], unname(p.adjust(expected_rr, "BH")[rownames(all)]))
})

test_that("a missing relative-risk p-value stays missing on its own row", {
  # A relative risk of 0 or below (or Inf) has no log-scale p-value. That row
  # must show NA in the relative-risk p-value columns, and the other rows'
  # adjusted values are computed as if it were still one of the tests.
  vims = list(a = make_vim("a", theta = 0.3, var_ic = 1, theta_rr = 1.1, var_ic_log_rr = 1),
              b = make_vim("b", theta = 0.2, var_ic = 1, theta_rr = Inf, var_ic_log_rr = NA),
              c = make_vim("c", theta = 0.1, var_ic = 1, theta_rr = 1.5, var_ic_log_rr = 1))
  out = varimpact:::compile_results(colnames_numeric = character(0),
                                    colnames_factor = names(vims),
                                    vim_numeric = list(),
                                    vim_factor = unname(vims),
                                    V = 2L)
  raw = out$results_raw
  expect_true(all(is.na(raw["b", c("rr_rawp", "rr_Holm", "rr_BH")])))
  expect_false(anyNA(raw[c("a", "c"), c("rr_rawp", "rr_Holm", "rr_BH")]))
  expected_rr = c(a = p_rr(1.1, 1), c = p_rr(1.5, 1))
  expect_equal(raw[c("a", "c"), "rr_rawp"], unname(expected_rr))
  expect_equal(raw[c("a", "c"), "rr_Holm"], unname(p.adjust(expected_rr, "holm", n = 3)))
  expect_equal(raw[c("a", "c"), "rr_BH"], unname(p.adjust(expected_rr, "BH", n = 3)))
})
