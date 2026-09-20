library(varimpact)
library(testthat)

context("adjust_pvalues")

# adjust_pvalues() replaced multtest::mt.rawp2adjp(p, c("Holm", "BH")). This
# is that function's Holm and BH arithmetic, transcribed from multtest so the
# equivalence is checked here without the Bioconductor dependency. It returns
# the adjusted values in the order of p, where multtest returned them sorted
# together with the sort index.
multtest_reference = function(rawp) {
  m = length(rawp)
  index = order(rawp)
  spval = rawp[index]

  holm = spval
  holm[1] = min(m * spval[1], 1)
  for (i in 2:m) {
    holm[i] = max(holm[i - 1], min((m - i + 1) * spval[i], 1))
  }

  bh = spval
  for (i in (m - 1):1) {
    bh[i] = min(bh[i + 1], min((m / i) * spval[i], 1, na.rm = TRUE), na.rm = TRUE)
    if (is.na(spval[i])) bh[i] = NA
  }

  out = cbind(rawp = spval, Holm = holm, BH = bh)
  out[order(index), , drop = FALSE]
}

test_that("Holm and BH match the multtest arithmetic", {
  set.seed(1, "L'Ecuyer-CMRG")
  for (m in c(2L, 3L, 7L, 40L)) {
    p = runif(m)^3
    expect_equal(varimpact:::adjust_pvalues(p), multtest_reference(p),
                 tolerance = 0, check.attributes = FALSE, info = m)
  }
})

test_that("known values", {
  p = c(0.01, 0.04, 0.03, 0.20)
  adjusted = varimpact:::adjust_pvalues(p)
  expect_identical(colnames(adjusted), c("rawp", "Holm", "BH"))
  expect_equal(adjusted[, "rawp"], p)
  # Holm: 4 * 0.01, then max(0.04, 3 * 0.03), max(0.09, 2 * 0.04), max(0.09, 0.2).
  expect_equal(adjusted[, "Holm"], c(0.04, 0.09, 0.09, 0.20))
  # BH, from the largest down: 0.2, then min(0.2, 4/3 * 0.04),
  # min(that, 4/2 * 0.03), min(that, 4/1 * 0.01).
  expect_equal(adjusted[, "BH"], c(0.04, 4 / 3 * 0.04, 4 / 3 * 0.04, 0.20))
})

test_that("adjustment never decreases a p-value and never exceeds one", {
  set.seed(2, "L'Ecuyer-CMRG")
  p = c(runif(20), 1, 0)
  adjusted = varimpact:::adjust_pvalues(p)
  expect_true(all(adjusted[, "Holm"] >= p))
  expect_true(all(adjusted[, "BH"] >= p))
  expect_true(all(adjusted[, c("Holm", "BH")] <= 1))
  # Holm is at least as conservative as BH.
  expect_true(all(adjusted[, "Holm"] >= adjusted[, "BH"]))
})

test_that("a missing p-value stays missing and still counts as a test", {
  p = c(0.01, NA, 0.04, NaN)
  adjusted = varimpact:::adjust_pvalues(p)
  expect_true(all(is.na(adjusted[c(2, 4), c("Holm", "BH")])))
  # The multipliers use all four p-values, as multtest did, not the two
  # that are observed. The reference handles NA the same way.
  expect_equal(adjusted[c(1, 3), "Holm"], c(0.04, 0.12))
  expect_equal(adjusted[c(1, 3), "BH"], c(0.04, 0.08))
  expect_equal(adjusted[, c("Holm", "BH")],
               multtest_reference(p)[, c("Holm", "BH")],
               check.attributes = FALSE)
})
