library(varimpact)
library(testthat)

context("cluster-robust variance")

# Shared fixture: families nested in sites, plus a crossed school variable.
set.seed(11, "L'Ecuyer-CMRG")
n_fam = 150
fam_size = sample(1:3, n_fam, replace = TRUE)
family = rep(seq_len(n_fam), fam_size)
n = length(family)
site = rep(sample(1:15, n_fam, replace = TRUE), fam_size)
school = sample(1:25, n, replace = TRUE)
# Influence-curve-like values with a family effect and a site effect.
ic = rnorm(n_fam)[family] + rnorm(15, 0, 0.5)[site] + rnorm(n)
ic = ic - mean(ic)

cluster_var = varimpact:::cluster_var
check_cluster_id = varimpact:::check_cluster_id

test_that("no clustering gives the ordinary variance exactly", {
  expect_identical(cluster_var(ic), var(ic))
  expect_identical(cluster_var(ic, NULL), var(ic))
  # One observation per cluster is the same as no clustering.
  expect_equal(cluster_var(ic, seq_len(n)), var(ic))
  expect_equal(cluster_var(ic, sample(seq_len(n))), var(ic))
})

test_that("one-way clustering matches the formula written out", {
  G = length(unique(site))
  s = tapply(ic - mean(ic), site, sum)
  expected = G / (G - 1) * sum(s^2) / n
  expect_equal(cluster_var(ic, site), expected)
  # Labels do not matter, only the partition.
  expect_equal(cluster_var(ic, paste0("s", site)), expected)
  expect_equal(cluster_var(ic, factor(site, levels = rev(unique(site)))),
               expected)
})

test_that("positive within-cluster correlation increases the variance", {
  expect_gt(cluster_var(ic, family), var(ic))
  expect_gt(cluster_var(ic, site), cluster_var(ic, family))
})

test_that("one-way and multiway variances match sandwich::vcovCL", {
  skip_if_not_installed("sandwich")
  # For an intercept-only model the coefficient's influence curve is y - ybar,
  # so vcovCL's clustered variance is cluster_var(ic, id) / n.
  y = ic
  fit = stats::lm(y ~ 1)
  vcl = function(...) {
    c(sandwich::vcovCL(fit, type = "HC0", cadjust = TRUE, ...))
  }
  expect_equal(cluster_var(ic, site) / n, vcl(cluster = ~site))
  expect_equal(cluster_var(ic, family) / n, vcl(cluster = ~family))
  expect_equal(cluster_var(ic, data.frame(site, school)) / n,
               vcl(cluster = ~site + school, multi0 = FALSE))
  # Three-way inclusion-exclusion.
  grade = sample(1:6, n, replace = TRUE)
  expect_equal(cluster_var(ic, data.frame(site, school, grade)) / n,
               vcl(cluster = ~site + school + grade, multi0 = FALSE))
})

test_that("a nested variable reduces to clustering on the outer variable", {
  expect_equal(cluster_var(ic, data.frame(family, site)),
               cluster_var(ic, site))
  expect_equal(cluster_var(ic, data.frame(site, family)),
               cluster_var(ic, site))
})

test_that("a non-positive multiway variance falls back to the largest one-way", {
  # A small case where V_A + V_B - V_AB < 0.
  ic_small = c(-0.9, 0.2, 1.6, -1.1, -0.1, 0.1)
  a = c(1, 1, 2, 2, 3, 3)
  b = c(1, 2, 1, 2, 1, 2)
  expect_warning(v <- cluster_var(ic_small, data.frame(a, b)),
                 "not positive")
  one_way = varimpact:::one_way_var
  expect_equal(v, max(one_way(ic_small, a), one_way(ic_small, b)))
})

test_that("check_cluster_id accepts vectors, data frames, matrices and lists", {
  v = check_cluster_id(site, n)
  expect_s3_class(v, "data.frame")
  expect_identical(names(v), "id")
  expect_true(is.integer(v$id))

  d = check_cluster_id(data.frame(site = site, school = school), n)
  expect_identical(names(d), c("site", "school"))

  m = check_cluster_id(cbind(site, school), n)
  expect_identical(names(m), c("site", "school"))

  l = check_cluster_id(list(site = site, school = school), n)
  expect_identical(names(l), c("site", "school"))

  # Unnamed columns still get distinct names.
  u = check_cluster_id(unname(cbind(site, school)), n)
  expect_equal(ncol(u), 2L)
  expect_false(anyDuplicated(names(u)) > 0)

  # Character ids are fine.
  expect_equal(cluster_var(ic, as.character(site)), cluster_var(ic, site))
})

test_that("check_cluster_id drops nested and singleton columns", {
  expect_message(d <- check_cluster_id(data.frame(family, site), n,
                                       verbose = TRUE),
                 "'family' is nested in 'site'")
  expect_identical(names(d), "site")

  expect_message(d <- check_cluster_id(data.frame(row = seq_len(n), site),
                                       n, verbose = TRUE),
                 "one observation per cluster")
  expect_identical(names(d), "site")

  # Only singletons: no clustering at all.
  expect_null(check_cluster_id(seq_len(n), n))
  expect_null(check_cluster_id(NULL, n))

  # Duplicate columns keep one copy.
  d = check_cluster_id(data.frame(a = site, b = site), n)
  expect_equal(ncol(d), 1L)

  # Crossed variables are both kept.
  d = check_cluster_id(data.frame(site, school), n)
  expect_identical(names(d), c("site", "school"))
})

test_that("check_cluster_id rejects bad input", {
  expect_error(check_cluster_id(site[-1], n), "rows but there are")
  bad = site
  bad[3] = NA
  expect_error(check_cluster_id(bad, n), "missing values")
  expect_error(check_cluster_id(data.frame(site, school = c(NA, school[-1])), n),
               "missing values")
})

test_that("ic_var matches cluster ids by row, not by position", {
  ic_var = varimpact:::ic_var
  cid = check_cluster_id(data.frame(site, school), n)
  rows = sample(seq_len(n), 200)
  sub_ic = ic[rows]
  expect_equal(ic_var(sub_ic, rows, cid),
               cluster_var(sub_ic, cid[rows, , drop = FALSE]))
  # Reordering the fold's observations does not change the variance.
  perm = sample(seq_along(rows))
  expect_equal(ic_var(sub_ic[perm], rows[perm], cid),
               ic_var(sub_ic, rows, cid))
})

test_that("ic_var keeps the unclustered path and guards its inputs", {
  ic_var = varimpact:::ic_var
  cid = check_cluster_id(site, n)
  expect_identical(ic_var(ic[1:50]), var(ic[1:50]))
  expect_identical(ic_var(ic[1:50], 1:50, NULL), var(ic[1:50]))
  expect_identical(ic_var(1, 1L, cid), NA_real_)
  expect_identical(ic_var(numeric(0), integer(0), cid), NA_real_)
  expect_error(ic_var(ic[1:10], 1:9, cid), "cannot be matched")
  expect_error(ic_var(ic[1:10], NULL, cid), "cannot be matched")
  # Every row from one cluster: too few clusters to estimate a variance.
  one_site = which(site == site[1])
  skip_if(length(one_site) < 2L)
  expect_warning(v <- ic_var(ic[one_site], one_site, cid), "fewer than 2")
  expect_identical(v, NA_real_)
})

test_that("cluster_df uses the variable with the fewest clusters", {
  cluster_df = varimpact:::cluster_df
  expect_identical(cluster_df(NULL), Inf)
  expect_equal(cluster_df(check_cluster_id(site, n)),
               length(unique(site)) - 1)
  expect_equal(cluster_df(check_cluster_id(data.frame(site, school), n)),
               min(length(unique(site)), length(unique(school))) - 1)
  # Nested family is dropped, so df comes from site.
  expect_equal(cluster_df(check_cluster_id(data.frame(family, site), n)),
               length(unique(site)) - 1)
})

test_that("cluster_sl_folds caps SuperLearner folds at the number of ids", {
  f = varimpact:::cluster_sl_folds
  expect_identical(f(1:100, 10), 10L)
  expect_identical(f(rep(1:4, 25), 10), 4L)
  expect_identical(f(rep(1, 10), 10), 1L)
})
