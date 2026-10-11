library(varimpact)
library(testthat)

context("create_cv_folds with clusters")

check_cluster_id = varimpact:::check_cluster_id
cluster_fold_groups = varimpact:::cluster_fold_groups
create_cv_folds = varimpact:::create_cv_folds

set.seed(21, "L'Ecuyer-CMRG")
n_fam = 120
fam_size = sample(1:4, n_fam, replace = TRUE)
family = rep(seq_len(n_fam), fam_size)
n = length(family)
site = rep(sample(1:12, n_fam, replace = TRUE), fam_size)
school = sample(1:20, n, replace = TRUE)
Y = rbinom(n, 1, 0.3)

intact = function(folds, g) {
  all(tapply(folds, g, function(x) length(unique(x))) == 1L)
}

test_that("without an id, fold assignment is unchanged", {
  Y2 = c(rep(0, 7), rep(1, 5))
  expect_identical(create_cv_folds(V = 3L, Y = Y2, id = NULL),
                   create_cv_folds(V = 3L, Y = Y2))
  expect_null(cluster_fold_groups(NULL, 5L))
})

test_that("every cluster falls within a single fold", {
  for (V in c(2L, 5L)) {
    folds = create_cv_folds(V, Y, id = family)
    expect_length(folds, n)
    expect_false(anyNA(folds))
    expect_setequal(unique(folds), seq_len(V))
    expect_true(intact(folds, family))
  }
})

test_that("fold sizes and events are close to balanced", {
  V = 5L
  folds = create_cv_folds(V, Y, id = family)
  sizes = tabulate(folds, nbins = V)
  # Greedy placement leaves folds within one largest cluster of each other.
  expect_lte(max(sizes) - min(sizes), max(fam_size))
  events = tabulate(folds[Y == 1], nbins = V)
  expect_lte(max(events) - min(events), max(tapply(Y, family, sum)))
})

test_that("a continuous outcome is also grouped", {
  Yc = rnorm(n)
  folds = create_cv_folds(4L, Yc, id = family)
  expect_true(intact(folds, family))
  expect_setequal(unique(folds), 1:4)
})

test_that("missing outcomes are placed with their cluster", {
  Ym = Y
  Ym[sample(n, 20)] = NA
  folds = create_cv_folds(3L, Ym, id = family)
  expect_false(anyNA(folds))
  expect_true(intact(folds, family))
})

test_that("grouped folds report errors clearly", {
  expect_error(create_cv_folds(5L, Y[1:10], id = rep(1:2, 5)),
               "fewer than V")
  expect_error(create_cv_folds(2L, Y, id = family[-1]),
               "cluster grouping has")
})

test_that("verbose prints the grouped fold breakdown", {
  expect_output(create_cv_folds(2L, Y, verbose = TRUE, id = family),
                "Cluster-grouped cross-validation")
})

test_that("nested ids are folded on the outer variable", {
  cid = check_cluster_id(data.frame(family, site), n)
  groups = cluster_fold_groups(cid, 5L)
  folds = create_cv_folds(5L, Y, id = groups)
  expect_true(intact(folds, site))
  expect_true(intact(folds, family))
})

test_that("crossed ids that chain together fall back with a warning", {
  cid = check_cluster_id(data.frame(site, school), n)
  expect_warning(groups <- cluster_fold_groups(cid, 5L), "crossed")
  # The coarsest variable with at least V clusters is site (12 < 20).
  expect_identical(groups, cid$site)
})

test_that("folds_by chooses the variable to fold on", {
  cid = check_cluster_id(data.frame(site, school), n)
  expect_identical(cluster_fold_groups(cid, 5L, folds_by = "school"),
                   cid$school)
  expect_error(cluster_fold_groups(cid, 5L, folds_by = "family"),
               "not one of the cluster variables")
  expect_error(cluster_fold_groups(NULL, 5L, folds_by = "site"),
               "no cluster id")
})

test_that("connected components keep separable crossed ids intact", {
  # Two blocks that never share a cluster on either variable, each repeated,
  # so that all ids can be kept whole.
  a = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6)
  b = c(1, 2, 1, 2, 3, 4, 3, 4, 5, 6, 5, 6)
  cid = check_cluster_id(data.frame(a, b), length(a))
  groups = cluster_fold_groups(cid, 3L)
  expect_identical(length(unique(groups)), 3L)
  folds = create_cv_folds(3L, rep(0:1, 6), id = groups)
  expect_true(intact(folds, a))
  expect_true(intact(folds, b))
})

test_that("an error is raised when no variable has V clusters", {
  cid = check_cluster_id(data.frame(a = rep(1:2, 10), b = rep(1:3, length.out = 20)), 20)
  expect_error(suppressWarnings(cluster_fold_groups(cid, 5L)),
               "No cluster variable has at least")
})
