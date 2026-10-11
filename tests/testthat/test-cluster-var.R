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

test_that("the leverage-corrected (CR3) variance matches vcovCL's HC3", {
  skip_if_not_installed("sandwich")
  # Unequal clusters: site 1 holds about a third of the rows.
  big = c(rep(1L, 150), site[-(1:150)] + 1L)
  y = ic
  fit = stats::lm(y ~ 1)
  w = rep(1, n)
  hc3 = function(...) {
    c(sandwich::vcovCL(fit, type = "HC3", cadjust = FALSE, ...))
  }
  expect_equal(cluster_var(ic, big, w) / n, hc3(cluster = ~big))
  expect_equal(cluster_var(ic, data.frame(big, school), w) / n,
               hc3(cluster = ~big + school, multi0 = FALSE))
  grade = sample(1:6, n, replace = TRUE)
  expect_equal(cluster_var(ic, data.frame(big, school, grade), w) / n,
               hc3(cluster = ~big + school + grade, multi0 = FALSE))
})

test_that("CR3 equals CR1 for equal clusters and the jackknife for a mean", {
  one_way = varimpact:::one_way_var
  g = rep(1:10, each = 6)
  x = ic[1:60]
  expect_equal(one_way(x, g, rep(1, 60)), one_way(x, g))
  # Unequal weights: the cluster jackknife of a weighted mean, centered at
  # the full-sample estimate.
  g = c(rep(1L, 30), rep(2:12, each = 4))
  m = length(g)
  w = seq_len(m) / sum(seq_len(m))
  y = ic[seq_len(m)]
  psi = sum(w * y)
  jack = vapply(1:12, function(k) sum((w * y)[g != k]) / sum(w[g != k]),
                numeric(1))
  expect_equal(one_way(m * w * (y - psi), g, w) / m,
               (11 / 12) * sum((jack - psi)^2))
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

test_that("check_cluster_id handles empty and badly named input", {
  # No columns at all, e.g. id = character(0) naming no columns of data.
  expect_null(check_cluster_id(list(), n))
  expect_null(check_cluster_id(data.frame(row.names = seq_len(n)), n))
  # Duplicate names are replaced so each column stays addressable.
  d = check_cluster_id(stats::setNames(data.frame(site, school), c("g", "g")), n)
  expect_identical(names(d), c("id1", "id2"))
})

test_that("a clustered variance needs at least two clusters", {
  expect_error(cluster_var(ic, rep(1L, n)), "At least 2 clusters")
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
  # Every row from one cluster: too few clusters to estimate a variance. The
  # clustered standard error comes from cv_cluster_var(), so this is quiet.
  one_site = which(site == site[1])
  skip_if(length(one_site) < 2L)
  expect_warning(v <- ic_var(ic[one_site], one_site, cid), NA)
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

cv_cluster_var = varimpact:::cv_cluster_var

test_that("cv_cluster_var weights folds and clusters once over all of them", {
  cid = check_cluster_id(site, n)
  fold = ifelse(site <= 7, 1L, 2L)
  ic_list = split(ic, fold)
  rows_list = split(seq_len(n), fold)
  # Each fold's curve is scaled by n / (V n_v), then clustered as a whole.
  sizes = lengths(ic_list)
  scaled = numeric(n)
  weights = numeric(n)
  for (v in 1:2) {
    scaled[rows_list[[v]]] = ic_list[[v]] * n / (2 * sizes[v])
    weights[rows_list[[v]]] = 1 / (2 * sizes[v])
  }
  # Clustered with the leverage (CR3) correction, each row weighted by its
  # share of the average over folds.
  expect_equal(cv_cluster_var(ic_list, rows_list, cid),
               cluster_var(scaled, site, weights))
  # The order of the folds' rows does not matter, only their row numbers.
  perm = lapply(rows_list, function(r) sample(seq_along(r)))
  expect_equal(cv_cluster_var(Map(`[`, ic_list, perm), Map(`[`, rows_list, perm),
                              cid),
               cv_cluster_var(ic_list, rows_list, cid))
  # Too few observations in a fold, or rows that do not line up.
  expect_identical(cv_cluster_var(list(ic[1:10], ic[11]), list(1:10, 11L), cid),
                   NA_real_)
  expect_identical(cv_cluster_var(list(), list(), cid), NA_real_)
  expect_error(cv_cluster_var(list(ic[1:10], ic[11:20]), list(1:10, 11:19), cid),
               "cannot be matched")
})

test_that("cv_cluster_var counts a cluster that spans folds once", {
  # Folds on site, so schools are split across folds.
  cid = check_cluster_id(data.frame(site, school), n)
  fold = ifelse(site <= 7, 1L, 2L)
  expect_gt(max(tapply(fold, school, function(f) length(unique(f)))), 1L)
  ic_list = split(ic, fold)
  rows_list = split(seq_len(n), fold)
  sizes = lengths(ic_list)
  scaled = numeric(n)
  weights = numeric(n)
  for (v in 1:2) {
    scaled[rows_list[[v]]] = ic_list[[v]] * n / (2 * sizes[v])
    weights[rows_list[[v]]] = 1 / (2 * sizes[v])
  }
  expect_equal(cv_cluster_var(ic_list, rows_list, cid),
               cluster_var(scaled, data.frame(site, school), weights))
})

test_that("copying every row within its cluster leaves the standard error", {
  # Two copies of each row double n, so the variance must double for the
  # standard error, sqrt(variance / n), to stay the same (issue #18).
  fold = ifelse(site <= 7, 1L, 2L)
  copies = rep(seq_len(n), each = 2)
  for (id in list(site, data.frame(site, school))) {
    id_copies = if (is.data.frame(id)) id[copies, ] else id[copies]
    once = cv_cluster_var(split(ic, fold), split(seq_len(n), fold),
                          check_cluster_id(id, n))
    twice = cv_cluster_var(split(ic[copies], fold[copies]),
                           split(seq_len(2 * n), fold[copies]),
                           check_cluster_id(id_copies, 2 * n))
    expect_equal(twice, 2 * once)
  }
})

test_that("cv_cluster_vars keeps the variation of a fold with one cluster", {
  # One cluster alone in fold 1, the rest in fold 2. Centering within each
  # fold would leave fold 1 contributing nothing.
  g = c(rep(1L, 60), rep(2:21, each = 5))
  m = length(g)
  cid = check_cluster_id(g, m)
  fold = ifelse(g == 1L, 1L, 2L)
  rows = split(seq_len(m), fold)
  pooled = function(shift) {
    curves = lapply(rows, function(r) {
      x = ic[r] + shift
      x - mean(x)
    })
    list(influence_curves = unname(curves),
         thetas = vapply(rows, function(r) 0.5 + shift + mean(ic[r]), numeric(1)),
         rows = unname(rows))
  }
  lo = pooled(0)
  hi = pooled(0.2)
  hi$influence_curves[[1]] = hi$influence_curves[[1]] + 0.1
  v = varimpact:::cv_cluster_vars(lo, hi, 2L, cid)
  expect_true(is.finite(v$rd) && v$rd > 0)
  expect_true(is.finite(v$log_rr) && v$log_rr > 0)
  # The risk difference uses the curves re-centered at the average estimate.
  recentered = lapply(1:2, function(k) {
    (hi$influence_curves[[k]] + hi$thetas[[k]] - mean(hi$thetas)) -
      (lo$influence_curves[[k]] + lo$thetas[[k]] - mean(lo$thetas))
  })
  expect_equal(v$rd, cv_cluster_var(recentered, lo$rows, cid))
})

test_that("grouped SuperLearner folds for g are stratified and keep clusters", {
  set.seed(12, "L'Ecuyer-CMRG")
  m = 300
  g = rep(1:60, each = 5)
  d = data.frame(A = rbinom(m, 1, 0.15), w = rnorm(m))
  fit = varimpact:::tmle_estimate_g(d, SL.library = c("SL.mean", "SL.glm"),
                                    id = g, V = 5)
  valid = fit$model$validRows
  expect_length(valid, 5L)
  fold = integer(m)
  for (k in seq_along(valid)) fold[valid[[k]]] = k
  expect_true(all(tapply(fold, g, function(x) length(unique(x))) == 1L))
  treated = tabulate(fold[d$A == 1], nbins = 5)
  expect_lte(max(treated) - min(treated), max(tapply(d$A, g, sum)))
  # Without stratification SuperLearner groups the folds itself.
  fit = varimpact:::tmle_estimate_g(d, SL.library = c("SL.mean", "SL.glm"),
                                    id = g, V = 5, stratify = FALSE)
  expect_identical(fit$type, "SuperLearner")
})

test_that("a single cluster falls back to ordinary SuperLearner folds", {
  set.seed(13, "L'Ecuyer-CMRG")
  m = 120
  W = data.frame(w = rnorm(m))
  A = rbinom(m, 1, 0.5)
  Y = rbinom(m, 1, plogis(W$w + A))
  fit = estimate_tmle2(Y = Y, A = A, W = W, family = "binomial",
                       Q.lib = c("SL.mean", "SL.glm"),
                       g.lib = c("SL.mean", "SL.glm"), V = 5,
                       Qbounds = c(0, 1), id = rep(1L, m))
  expect_length(fit$g_model$validRows, 2L)
  expect_length(fit$q_model$validRows, 2L)
})

test_that("factors that never cross cluster-grouped folds are set aside", {
  drop_cols = varimpact:::drop_cluster_level_columns
  d = data.frame(site_name = factor(paste0("s", site)),       # the id itself
                 member = factor(paste(family, seq_len(n))),  # finer than site
                 kind = factor(c("a", "b", "c")[site %% 3 + 1]),  # spans sites
                 size = (site * 7)[seq_len(n)],              # numeric, per site
                 x = ic)
  expect_identical(drop_cols(d, NULL), d)
  expect_message(kept <- drop_cols(d, site), "Not analyzing site_name, member")
  expect_identical(names(kept), c("kind", "size", "x"))
  # Character columns and matrices are handled the same way.
  m = cbind(site_name = paste0("s", site), kind = c("a", "b", "c")[site %% 3 + 1])
  expect_message(kept <- drop_cols(m, site), "Not analyzing site_name")
  expect_identical(colnames(kept), "kind")
  expect_error(suppressMessages(drop_cols(d[c("site_name", "member")], site)),
               "no variables left")
})

test_that("contrast_df counts the clusters behind the scarcer level", {
  cdf = varimpact:::contrast_df
  cid = check_cluster_id(data.frame(site, school), n)
  fold = ifelse(site <= 7, 1L, 2L)
  preds = function(rows, a, delta = rep(1, length(rows))) {
    list(val_preds = data.frame(A = a, delta = delta, row_index = rows))
  }
  count = function(rows) {
    c(length(unique(site[rows])), length(unique(school[rows])))
  }
  # Fold 1 picks sites 1-2 as its low level, fold 2 picks sites 14-15 as its
  # high level, and the other level is the rest of the fold.
  fold_results = list(
    list(level_min = preds(which(fold == 1), as.integer(site[fold == 1] <= 2)),
         level_max = preds(which(fold == 1), as.integer(site[fold == 1] > 2))),
    list(level_min = preds(which(fold == 2), as.integer(site[fold == 2] < 14)),
         level_max = preds(which(fold == 2), as.integer(site[fold == 2] >= 14))))
  expected = min(pmin(count(which(site <= 2)), count(which(site > 2 & fold == 1))) +
                   pmin(count(which(site < 14 & fold == 2)), count(which(site >= 14)))) - 1
  expect_equal(cdf(fold_results, cid), expected)
  expect_equal(cdf(fold_results, check_cluster_id(site, n)), 2 + 2 - 1)
  # Rows with delta = 0 do not count: site 2's rows become unobserved.
  p = fold_results[[1]]$level_min$val_preds
  p$delta[site[p$row_index] == 2] = 0
  fold_results[[1]]$level_min$val_preds = p
  expect_equal(cdf(fold_results, check_cluster_id(site, n)), 1 + 2 - 1)
  # A fold without predictions contributes nothing.
  fold_results[[2]] = list(level_min = list(val_preds = NULL),
                           level_max = list(val_preds = NULL))
  expect_equal(cdf(fold_results, check_cluster_id(site, n)), 1 - 1)
})
