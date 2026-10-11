library(varimpact)
library(testthat)

context("varimpact() with cluster-robust standard errors")

future::plan("sequential")

set.seed(41, "L'Ecuyer-CMRG")
n_fam = 70
fam_size = sample(1:3, n_fam, replace = TRUE)
family = rep(seq_len(n_fam), fam_size)
N = length(family)
site = rep(sample(1:10, n_fam, replace = TRUE), fam_size)
school = sample(1:12, N, replace = TRUE)
u = rnorm(n_fam, 0, 1)[family]
X = data.frame(x1 = rnorm(n_fam)[family] + rnorm(N, 0, 0.3),
               x2 = rnorm(N),
               f1 = factor(sample(c("a", "b", "c"), N, replace = TRUE)))
Y = rbinom(N, 1, plogis(X$x1 + u))

run = function(...) {
  suppressWarnings(
    varimpact(Y = Y, data = X, V = 2L,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"),
              bins_numeric = 3L, ...))
}

test_that("an id with one observation per cluster changes nothing", {
  set.seed(1, "L'Ecuyer-CMRG")
  base = run()
  set.seed(1, "L'Ecuyer-CMRG")
  expect_message(singleton <- run(id = seq_len(N)), "one observation per cluster")
  expect_identical(singleton$results_all, base$results_all)
  expect_identical(singleton$cv_folds, base$cv_folds)
  expect_null(singleton$cluster_id)
  expect_identical(singleton$cluster_df, Inf)
  expect_false(any(grepl("Cluster-robust", capture.output(print(base)))))
})

test_that("clustering by family keeps families in one fold and uses t df", {
  set.seed(2, "L'Ecuyer-CMRG")
  vim = run(id = family)
  expect_true(all(tapply(vim$cv_folds, family,
                         function(x) length(unique(x))) == 1L))
  expect_equal(vim$cluster_df, n_fam - 1)
  expect_identical(names(vim$cluster_id), "id")
  expect_false(is.null(vim$results_all))
  expect_true(all(is.finite(vim$results_all$Estimate)))
  expect_output(print(vim), paste0("Cluster-robust standard errors: id \\(",
                                   n_fam, " clusters\\)\\. Intervals and ",
                                   "p-values use a t distribution with the df ",
                                   "shown for each variable \\(at most ",
                                   n_fam - 1, "\\)"))
  # Each variable's df sits next to its interval, capped at the global df.
  expect_identical(names(vim$results_all)[1:4],
                   c("Type", "Estimate", "CI95", "df"))
  expect_true(all(vim$results_all$df <= n_fam - 1))
  expect_identical(vim$results_raw$df, vim$results_all$df)
})

test_that("nested family and site reduce to site", {
  set.seed(3, "L'Ecuyer-CMRG")
  expect_message(vim <- run(id = data.frame(family = family, site = site)),
                 "nested in 'site'")
  expect_identical(names(vim$cluster_id), "site")
  expect_equal(vim$cluster_df, length(unique(site)) - 1)
  expect_true(all(tapply(vim$cv_folds, site,
                         function(x) length(unique(x))) == 1L))
})

test_that("crossed ids run with multiway variances and folds_by", {
  set.seed(4, "L'Ecuyer-CMRG")
  vim = run(id = data.frame(site = site, school = school), folds_by = "site")
  expect_identical(names(vim$cluster_id), c("site", "school"))
  expect_equal(vim$cluster_df,
               min(length(unique(site)), length(unique(school))) - 1)
  expect_true(all(tapply(vim$cv_folds, site,
                         function(x) length(unique(x))) == 1L))
})

test_that("a cluster holding most of the data still yields standard errors", {
  # One family holds over half of the rows. It shares its fold with another
  # family rather than filling the fold alone, and every variable gets an
  # interval and a p-value.
  big = family
  big[family <= 40] = 1L
  set.seed(8, "L'Ecuyer-CMRG")
  vim = run(id = big)
  per_fold = tapply(big, vim$cv_folds, function(x) length(unique(x)))
  expect_true(all(per_fold >= 2L))
  expect_identical(nrow(vim$results_all), 3L)
  expect_true(all(grepl("^\\(", vim$results_all$CI95)))
  expect_false(anyNA(vim$results_all$`P-value`))
})

test_that("a family-level factor in few families gets fewer df", {
  # A flag set for whole families, in about one family in four: its contrast
  # rests on those families, not on all of them.
  flagged = unique(family)[seq(1, n_fam, by = 4)]
  Xf = cbind(X, flag = factor(family %in% flagged))
  set.seed(10, "L'Ecuyer-CMRG")
  vim = suppressWarnings(
    varimpact(Y = Y, data = Xf, id = family, V = 2L,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L))
  expect_true("flag" %in% rownames(vim$results_all))
  expect_lte(vim$results_all["flag", "df"], length(flagged) - 1)
  expect_gt(vim$results_all["x1", "df"], vim$results_all["flag", "df"])
})

test_that("a continuous outcome with id gives results on the outcome scale", {
  # Rescaled to the outcome's range internally; the clustered variance must be
  # computed on the original scale, like the estimate.
  Yc = 10 + 3 * X$x1 + rnorm(n_fam, 0, 2)[family] + rnorm(N)
  set.seed(11, "L'Ecuyer-CMRG")
  vim = suppressWarnings(
    varimpact(Y = Yc, data = X, id = family, family = "gaussian", V = 2L,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L))
  res = vim$results_all
  expect_true(all(c("x1", "x2", "f1") %in% rownames(res)))
  # x1 moves the outcome by several units between its low and high bins.
  expect_gt(res["x1", "Estimate"], 1)
  expect_lt(res["x1", "P-value"], 0.05)
  expect_true(all(is.finite(res$df)) && all(res$df <= n_fam - 1))
  bounds = regmatches(res["x1", "CI95"], gregexpr("-?[0-9.]+", res["x1", "CI95"]))
  bounds = as.numeric(bounds[[1]])
  expect_true(bounds[1] < res["x1", "Estimate"] && res["x1", "Estimate"] < bounds[2])
})

test_that("missing outcomes are handled with id", {
  # Rows missing Y stay in their cluster's fold but add nothing to the
  # residual term, and do not count toward a level's clusters.
  Ym = Y
  Ym[seq(1, N, by = 9)] = NA
  set.seed(12, "L'Ecuyer-CMRG")
  vim = suppressWarnings(
    varimpact(Y = Ym, data = X, id = family, V = 2L,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L))
  res = vim$results_all
  expect_true(all(c("x1", "x2", "f1") %in% rownames(res)))
  expect_false(anyNA(res$`P-value`))
  expect_true(all(is.finite(res$df)) && all(res$df <= n_fam - 1))
  expect_true(all(tapply(vim$cv_folds, family,
                         function(x) length(unique(x))) == 1L))
})

test_that("duplicating every row under a shared id leaves the intervals", {
  # Issue #18: two copies of each row, sharing an id, carry no more
  # information than one. Without id the intervals narrow by about sqrt(2);
  # with id they should match the original data's. Point estimates can still
  # move, since the folds, and so each variable's low and high levels, differ.
  set.seed(18, "L'Ecuyer-CMRG")
  n1 = 200
  Xi = data.frame(x1 = rnorm(n1), x2 = rnorm(n1),
                  f1 = factor(sample(c("a", "b", "c"), n1, replace = TRUE)))
  Yi = rbinom(n1, 1, plogis(-0.5 + 0.8 * Xi$x1 + 0.5 * (Xi$f1 == "b")))
  copies = rep(seq_len(n1), each = 2)
  half_width = function(...) {
    set.seed(19, "L'Ecuyer-CMRG")
    res = suppressWarnings(suppressMessages(
      varimpact(V = 2L, Q.library = c("SL.mean", "SL.glm"),
                g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L, ...)
    ))$results_all
    bounds = strsplit(gsub("[()]", "", res$CI95), " - ", fixed = TRUE)
    stats::setNames(vapply(bounds, function(b) diff(as.numeric(b)) / 2,
                           numeric(1)),
                    rownames(res))
  }
  vars = c("x1", "x2", "f1")
  original = half_width(Y = Yi, data = Xi)[vars]
  copied = half_width(Y = Yi[copies], data = Xi[copies, ])[vars]
  copied_id = half_width(Y = Yi[copies], data = Xi[copies, ], id = copies)[vars]
  expect_false(anyNA(c(original, copied, copied_id)))
  # Across data seeds the ratios ran 0.92-1.21 with id and 1.39-1.8 without.
  expect_true(all(copied_id / original > 1 / 1.3 & copied_id / original < 1.3))
  expect_true(all(copied_id / copied > 1.3))
})

test_that("a cluster id left in data as a factor is set aside with a message", {
  # site is the cluster and is also a column of data. Each site sits in one
  # fold, so its levels never match between training and validation.
  Xs = cbind(X, site_f = factor(site))
  set.seed(9, "L'Ecuyer-CMRG")
  expect_message(
    vim <- suppressWarnings(
      varimpact(Y = Y, data = Xs, id = site, V = 2L,
                Q.library = c("SL.mean", "SL.glm"),
                g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L)),
    "Not analyzing site_f")
  expect_false("site_f" %in% rownames(vim$results_all))
  expect_true(all(c("x1", "x2", "f1") %in% rownames(vim$results_all)))
})

test_that("bad cluster input is rejected before any fitting", {
  expect_error(varimpact(Y = Y, data = X, id = family[-1]),
               "rows but there are")
  bad = family
  bad[1] = NA
  expect_error(varimpact(Y = Y, data = X, id = bad), "missing values")
  expect_error(varimpact(Y = Y, data = X, id = family, folds_by = "site"),
               "not one of the cluster variables")
  expect_error(varimpact(Y = Y[-1], data = X), "Y has")
})

test_that("id can name columns of data, which are then not analyzed", {
  Xid = cbind(X, family = family, site = site)

  set.seed(5, "L'Ecuyer-CMRG")
  by_name = suppressWarnings(
    varimpact(Y = Y, data = cbind(X, family = family), id = "family", V = 2L,
              Q.library = c("SL.mean", "SL.glm"),
              g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L))
  set.seed(5, "L'Ecuyer-CMRG")
  by_vector = run(id = family)
  # Same result as passing the vector itself, and family is not a predictor.
  expect_identical(by_name$results_all, by_vector$results_all)
  expect_false("family" %in% rownames(by_name$results_all))

  # Several names: multiway (here nested, so site is kept).
  set.seed(6, "L'Ecuyer-CMRG")
  expect_message(
    two <- suppressWarnings(
      varimpact(Y = Y, data = Xid, id = c("family", "site"), V = 2L,
                Q.library = c("SL.mean", "SL.glm"),
                g.library = c("SL.mean", "SL.glm"), bins_numeric = 3L)),
    "nested in 'site'")
  expect_identical(names(two$cluster_id), "site")
  expect_false(any(c("family", "site") %in% rownames(two$results_all)))
})

test_that("a character id as long as the data is used as ids, not names", {
  set.seed(7, "L'Ecuyer-CMRG")
  a = run(id = paste0("fam", family))
  set.seed(7, "L'Ecuyer-CMRG")
  b = run(id = family)
  # Equal rather than identical: relabeling the clusters changes the order in
  # which their sums are added, which moves the last bits on platforms without
  # extended-precision accumulation (Apple silicon).
  expect_equal(a$results_all, b$results_all)
})

test_that("an id that is NULL by mistake warns that clustering is off", {
  Xf = cbind(X, family = family)
  # family = "poisson" stops right after the id checks, before any fitting.
  expect_warning(
    expect_error(varimpact(Y = Y, data = X, id = Xf$famly, family = "poisson"),
                 "Family must be"),
    "NOT clustered")
  expect_warning(
    expect_error(varimpact(Y = Y, data = X, id = Xf[["famly"]],
                           family = "poisson"),
                 "Family must be"),
    "NOT clustered")
  # Leaving id out, passing NULL on purpose, or a wrapper forwarding its own
  # NULL default does not warn.
  wrapper = function(Y, data, id = NULL) {
    varimpact(Y = Y, data = data, id = id, family = "poisson")
  }
  for (call in list(quote(varimpact(Y = Y, data = X, family = "poisson")),
                    quote(varimpact(Y = Y, data = X, id = NULL,
                                    family = "poisson")),
                    quote(wrapper(Y, X)))) {
    expect_warning(expect_error(eval(call), "Family must be"), NA)
  }
})

test_that("verbose output says which clusters are used", {
  # family = "poisson" stops right after the cluster setup, before any fitting.
  expect_output(
    expect_error(varimpact(Y = Y, data = cbind(X, family = family),
                           id = "family", family = "poisson", verbose = TRUE),
                 "Family must be"),
    paste0("Clustering on column\\(s\\) family of data.*",
           "Cluster-robust standard errors using: family \\(", n_fam,
           " clusters\\).*t critical values with at most ", n_fam - 1))
})

test_that("unknown id column names are reported clearly", {
  expect_error(varimpact(Y = Y, data = X, id = "familly"),
               "not found in data: familly")
  expect_error(varimpact(Y = Y, data = X["x1"], id = "x1"),
               "no variables left")
})
