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
  expect_identical(a$results_all, b$results_all)
})

test_that("an id that is NULL by mistake warns that clustering is off", {
  Xf = cbind(X, family = family)
  # family = "poisson" stops right after the id checks, before any fitting.
  expect_warning(
    expect_error(varimpact(Y = Y, data = X, id = Xf$famly, family = "poisson"),
                 "Family must be"),
    "NOT clustered")
  # Leaving id out entirely does not warn.
  expect_warning(
    expect_error(varimpact(Y = Y, data = X, family = "poisson"),
                 "Family must be"),
    NA)
})

test_that("unknown id column names are reported clearly", {
  expect_error(varimpact(Y = Y, data = X, id = "familly"),
               "not found in data: familly")
  expect_error(varimpact(Y = Y, data = X["x1"], id = "x1"),
               "no variables left")
})
