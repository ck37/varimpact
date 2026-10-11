library(varimpact)
library(testthat)

context("simulation study of cluster-robust inference (opt-in)")

# A simulation study of varimpact()'s cluster-robust standard errors. It is
# skipped unless VARIMPACT_RUN_SIMULATIONS=true, because it fits varimpact()
# a few thousand times:
#
#   VARIMPACT_RUN_SIMULATIONS=true VARIMPACT_SIM_CORES=4 \
#     Rscript -e 'devtools::test(filter = "cluster-simulation")'
#
# VARIMPACT_SIM_CORES runs replications in parallel with
# parallel::mclapply() (not on Windows). VARIMPACT_SIM_REPS sets the number
# of replications for every scenario; by default each scenario uses the
# number in its definition below, which reproduces the reference results.
# The full-pipeline test took about an hour on 3 cores; the test of the
# variance formula alone takes about a minute.
#
# Each replication draws clustered data with an outcome that has a strong
# cluster effect, and three predictors:
#   x_null  varies within clusters but also has a cluster component; no effect
#   f_null  a three-level factor constant within clusters; no effect
#   x_eff   has a cluster component and a real effect on the outcome
# and fits varimpact() with and without id (scenario F also with
# folds_by = "site"). A 95% interval that excludes 0
# (a two-sided p-value below 0.05) for x_null or f_null is a false positive,
# so those rates should be near 5%; with 300 replications, rates between
# about 2.5% and 7.5% are consistent with 5% (2% to 8% with 200).
#
# Reference results: see the end of this file.

simulate_clustered_data = function(G, size_range, dominant = 0, icc_scale = 1) {
  size = sample(size_range[1]:size_range[2], G, replace = TRUE)
  # Optionally make the first cluster hold a share `dominant` of the rows.
  if (dominant > 0) size[1] = round(dominant * sum(size[-1]) / (1 - dominant))
  g = rep(seq_len(G), size)
  N = length(g)
  u = rnorm(G, 0, icc_scale)[g]
  x_null = rnorm(G)[g] + rnorm(N, 0, 0.5)
  f_null = factor(sample(c("a", "b", "c"), G, replace = TRUE)[g])
  x_eff = rnorm(G, 0, 0.7)[g] + rnorm(N, 0, 0.7)
  Y = rbinom(N, 1, plogis(-0.3 + 0.6 * x_eff + u))
  list(Y = Y, X = data.frame(x_null, f_null, x_eff), g = g)
}

# Crossed clustering: every row belongs to one of 20 sites and, independently,
# one of 15 schools, and the outcome and the covariates vary with both. The
# two variables chain every row together, so no fold assignment can keep
# both intact: varimpact() folds on school (the coarser variable, with a
# warning), and each site spans both folds, so a training fold shares sites
# with its validation fold. f_null is constant within sites. This scenario is
# also fitted with folds_by = "site", which keeps sites intact instead.
simulate_crossed_data = function(n = 600, sites = 20, schools = 15) {
  site = sample(seq_len(sites), n, replace = TRUE)
  school = sample(seq_len(schools), n, replace = TRUE)
  u = rnorm(sites, 0, 0.7)[site] + rnorm(schools, 0, 0.7)[school]
  x_null = rnorm(sites)[site] + rnorm(schools)[school] + rnorm(n, 0, 0.5)
  f_null = factor(sample(c("a", "b", "c"), sites, replace = TRUE)[site])
  x_eff = rnorm(sites, 0, 0.5)[site] + rnorm(schools, 0, 0.5)[school] +
    rnorm(n, 0, 0.7)
  Y = rbinom(n, 1, plogis(-0.3 + 0.6 * x_eff + u))
  list(Y = Y, X = data.frame(x_null, f_null, x_eff),
       g = data.frame(site = site, school = school))
}

simulation_scenarios = list(
  list(id = 1, name = "A: 40 clusters", reps = 300,
       args = list(G = 40, size_range = c(5, 20))),
  list(id = 2, name = "B: 15 clusters", reps = 300,
       args = list(G = 15, size_range = c(10, 30))),
  list(id = 3, name = "C: one cluster holds 30%", reps = 200,
       args = list(G = 40, size_range = c(5, 20), dominant = 0.3)),
  list(id = 4, name = "D: no cluster effect on Y", reps = 200,
       args = list(G = 40, size_range = c(5, 20), icc_scale = 0)),
  list(id = 5, name = "E: one cluster holds 50%", reps = 200,
       args = list(G = 40, size_range = c(5, 20), dominant = 0.5)),
  list(id = 6, name = "F: crossed site x school", reps = 200,
       args = list(), simulate = simulate_crossed_data, folds_by = "site"))

simulation_replication = function(rep, scenario) {
  set.seed(1000 * scenario$id + rep, "L'Ecuyer-CMRG")
  simulate = if (is.null(scenario$simulate)) simulate_clustered_data else
    scenario$simulate
  d = do.call(simulate, scenario$args)
  fit = function(id, folds_by = NULL) {
    vim = tryCatch(suppressWarnings(suppressMessages({
      utils::capture.output(
        v <- varimpact(Y = d$Y, data = d$X, V = 2L, id = id,
                       folds_by = folds_by,
                       Q.library = c("SL.mean", "SL.glm"),
                       g.library = c("SL.mean", "SL.glm"),
                       bins_numeric = 4L, parallel = FALSE))
      v
    })), error = function(e) NULL)
    vars = c("x_null", "f_null", "x_eff")
    raw = vim$results_raw
    m = if (is.null(raw)) rep(NA_integer_, 3L) else match(vars, rownames(raw))
    data.frame(var = vars,
               estimate = if (is.null(raw)) NA_real_ else raw$AvePsi[m],
               p = if (is.null(raw)) NA_real_ else raw$rawp[m])
  }
  out = rbind(cbind(fit(NULL), analysis = "no id"),
              cbind(fit(d$g), analysis = "id"))
  # Optionally also fold on a chosen clustering variable.
  if (!is.null(scenario$folds_by)) {
    out = rbind(out, cbind(fit(d$g, scenario$folds_by),
                           analysis = paste("id, folds_by =", scenario$folds_by)))
  }
  out
}

run_cluster_simulation = function(scenarios, reps = NULL, cores = 1L) {
  if (.Platform$OS.type == "windows") cores = 1L
  do.call(rbind, lapply(scenarios, function(scenario) {
    n = if (is.null(reps)) scenario$reps else reps
    res = parallel::mclapply(seq_len(n), simulation_replication,
                             scenario = scenario, mc.cores = cores)
    res = do.call(rbind, res)
    res$scenario = scenario$name
    res
  }))
}

summarize_cluster_simulation = function(results) {
  results$reject = results$p < 0.025 | results$p > 0.975
  cells = split(results, list(results$scenario, results$var, results$analysis),
                drop = TRUE)
  out = do.call(rbind, lapply(cells, function(d) {
    ok = !is.na(d$p)
    data.frame(scenario = d$scenario[1], var = d$var[1],
               analysis = d$analysis[1], reps = nrow(d),
               estimated = mean(ok), reject = mean(d$reject[ok]))
  }))
  out[order(out$scenario, out$var, out$analysis), ]
}

test_that("with id, variables without an effect are flagged about 5% of the time", {
  skip_if_not(identical(Sys.getenv("VARIMPACT_RUN_SIMULATIONS"), "true"),
              "set VARIMPACT_RUN_SIMULATIONS=true to run the simulation study")
  future::plan("sequential")
  reps = Sys.getenv("VARIMPACT_SIM_REPS")
  reps = if (nzchar(reps)) as.integer(reps) else NULL
  cores = as.integer(Sys.getenv("VARIMPACT_SIM_CORES", "1"))

  summary = summarize_cluster_simulation(
    run_cluster_simulation(simulation_scenarios, reps, cores))
  print(summary, row.names = FALSE, digits = 3)

  # With clustering, each null variable's false-positive rate should be within
  # Monte Carlo error of 5%: at most 5% plus three standard errors. This holds
  # where every clustering variable stays within folds. With crossed
  # variables folded by default, a fold's training and validation data share
  # clusters, which biases variables that vary with them (see the reference
  # results at the end of this file); there id should still beat no id, and
  # folding on the variable they vary with should restore 5%.
  null = summary[summary$analysis != "no id" & summary$var != "x_eff", ]
  crossed_default = grepl("crossed", null$scenario) & null$analysis == "id"
  strict = null[!crossed_default, ]
  limit = 0.05 + 3 * sqrt(0.05 * 0.95 / strict$reps)
  expect_true(all(strict$reject <= limit),
              info = paste(strict$scenario, strict$analysis, strict$var,
                           round(strict$reject, 3), collapse = "; "))
  if (any(crossed_default)) {
    with_id = null[crossed_default, ]
    without = summary[summary$analysis == "no id" &
                        summary$scenario %in% with_id$scenario, ]
    without = without[match(paste(with_id$scenario, with_id$var),
                             paste(without$scenario, without$var)), ]
    expect_true(all(with_id$reject < without$reject))
  }
})

# A faster study of the variance formula alone, without fitting varimpact().
# The estimate mimics varimpact's: the average over folds of each fold's mean,
# with the influence curve each row's deviation. It compares three ways to
# get a clustered variance:
#   per_fold       the average of per-fold clustered variances
#   fold_centered  one clustered variance over all folds, with each fold's
#                  curve centered at that fold's estimate
#   varimpact      one clustered variance over all folds, centered at the
#                  overall estimate (cv_cluster_var(), as varimpact() does)
# and reports 95% interval coverage, which should be near 95%.
variance_design_coverage = function(make_y, fold_of, cluster_id, df, reps) {
  ic_var = varimpact:::ic_var
  cv_cluster_var = varimpact:::cv_cluster_var
  est = matrix(NA_real_, reps, 4L,
               dimnames = list(NULL, c("estimate", "per_fold",
                                       "fold_centered", "varimpact")))
  for (r in seq_len(reps)) {
    fold = fold_of()
    y = make_y()
    n = length(y)
    V = max(fold)
    means = tapply(y, fold, mean)
    rows = split(seq_len(n), fold)
    ic = y - means[fold]
    est[r, "estimate"] = mean(means)
    est[r, "per_fold"] = sqrt(mean(vapply(seq_len(V), function(v) {
      ic_var(ic[rows[[v]]], rows[[v]], cluster_id)
    }, numeric(1))) / n)
    est[r, "fold_centered"] = sqrt(cv_cluster_var(split(ic, fold), rows,
                                                  cluster_id) / n)
    est[r, "varimpact"] = sqrt(cv_cluster_var(split(y - mean(means), fold),
                                              rows, cluster_id) / n)
  }
  crit = stats::qt(0.975, df)
  covered = function(se) mean(abs(est[, "estimate"]) < crit * se, na.rm = TRUE)
  c(per_fold = covered(est[, "per_fold"]),
    fold_centered = covered(est[, "fold_centered"]),
    varimpact = covered(est[, "varimpact"]))
}

test_that("the clustered variance keeps 95% coverage across cluster designs", {
  skip_if_not(identical(Sys.getenv("VARIMPACT_RUN_SIMULATIONS"), "true"),
              "set VARIMPACT_RUN_SIMULATIONS=true to run the simulation study")
  reps = as.integer(Sys.getenv("VARIMPACT_SIM_REPS", "2000"))
  ccid = varimpact:::check_cluster_id
  folds_for = function(g, V) {
    function() varimpact:::assign_group_folds(g, V, rep(0, length(g)))
  }
  one_way = function(G, V, seed) {
    set.seed(seed)
    g = rep(seq_len(G), sample(5:40, G, replace = TRUE))
    list(make_y = function() rnorm(G, 0, 0.5)[g] + rnorm(length(g)),
         fold_of = folds_for(g, V), cluster_id = ccid(g, length(g)),
         df = G - 1)
  }
  designs = list(
    "crossed site x school, folds on site" = local({
      set.seed(5)
      n = 600
      site = sample(1:20, n, TRUE)
      school = sample(1:15, n, TRUE)
      list(make_y = function() {
             rnorm(20, 0, 0.3)[site] + rnorm(15, 0, 0.6)[school] + rnorm(n)
           },
           fold_of = function() ifelse(site <= 10, 1L, 2L),
           cluster_id = ccid(data.frame(site, school), n), df = 14)
    }),
    "one-way, 40 clusters" = one_way(40, 2L, 6),
    "one-way, 12 clusters" = one_way(12, 2L, 8),
    "one-way, 40 clusters, V = 5" = one_way(40, 5L, 9))
  # One cluster holding a growing share of the rows, beside 30 clusters of 8.
  for (share in c(0.05, 0.1, 0.2, 0.3, 0.5)) {
    big = round(share * 30 * 8 / (1 - share))
    g = c(rep(1L, big), rep(2:31, each = 8))
    designs[[sprintf("one cluster holds %d%%", round(100 * share))]] = local({
      g = g
      list(make_y = function() rnorm(31, 0, 0.5)[g] + rnorm(length(g)),
           fold_of = folds_for(g, 2L), cluster_id = ccid(g, length(g)),
           df = 30)
    })
  }

  set.seed(1)
  coverage = t(vapply(designs, function(d) {
    variance_design_coverage(d$make_y, d$fold_of, d$cluster_id, d$df, reps)
  }, numeric(3)))
  print(round(coverage, 3))

  # varimpact's variance stays near 95% in every design, including one
  # cluster holding half the data: at least 90%, less three Monte Carlo
  # standard errors.
  limit = 0.9 - 3 * sqrt(0.95 * 0.05 / reps)
  expect_true(all(coverage[, "varimpact"] >= limit),
              info = paste(rownames(coverage), round(coverage[, "varimpact"], 3),
                           collapse = "; "))
})

# Reference results, from running both tests above with VARIMPACT_SIM_CORES=3
# (October 2026, R 4.3, the code as of the commit that last updated them).
#
# Full pipeline: the share of replications in which each variable's 95%
# interval excluded 0, without and with id. For x_null and f_null that is
# the false-positive rate; for x_eff it is the rate of flagging a real effect
# (inflated without id, since that analysis also flags null variables).
#
#   scenario                      x_null         f_null         x_eff
#                                 no id  id      no id  id      no id  id
#   A: 40 clusters                19.3%  8.0%    42.0%  5.3%    90.3%  70.3%
#   B: 15 clusters                27.0%  5.1%    57.7%  0.0%    72.7%  36.0%
#   C: one cluster holds 30%      27.5%  9.0%    66.0%  3.5%    93.0%  65.0%
#   D: no cluster effect on Y      5.5%  7.0%    11.5%  6.5%    97.0%  96.0%
#   E: one cluster holds 50%      27.5%  3.5%    72.5%  1.6%    96.0%  53.8%
#   F: crossed site x school      22.5%  8.0%    36.5% 13.0%    96.0%  80.0%
#      with folds_by = "site"             7.0%           3.0%           82.5%
#
# With id, f_null was estimated in 98.3% of B's replications and 91.5% of
# E's (82% before every fold was given at least two clusters): when the
# large cluster's fold holds few clusters, the other fold's training data can
# hold a single level of f_null. With that rule, rates in A, B and D moved by
# at most half a point, and C's and E's by amounts within Monte Carlo error.
#
# In F, varimpact() folds on school by default, so each site spans both
# folds. With id, f_null's 13% is bias, not a too-small variance: its
# estimate averaged +0.069 (0.73 of its standard deviation) where the truth is
# 0, and nearly all rejections were in the upper tail. Each training fold
# picks f_null's most different levels partly because of chance site
# effects, and the same sites then confirm the difference in the validation
# fold. Folding on site keeps sites intact and removes that bias (f_null
# 3.0%); a variable that varied with school would then be exposed instead.
# varimpact() warns about this when it falls back to one variable.
#
# Variance formula alone (2,000 replications): 95% interval coverage.
#
#   design                                per_fold  fold_centered  varimpact
#   crossed site x school, folds on site  87.0%     94.3%          94.3%
#   one-way, 40 clusters                  94.5%     94.4%          94.6%
#   one-way, 12 clusters                  91.3%     90.9%          93.5%
#   one-way, 40 clusters, V = 5           93.7%     92.7%          94.5%
#   one cluster holds 5%                  94.6%     94.2%          94.5%
#   one cluster holds 10%                 93.5%     93.4%          94.3%
#   one cluster holds 20%                 89.3%     90.1%          93.1%
#   one cluster holds 30%                 80.2%     84.5%          93.4%
#   one cluster holds 50%                 37.4%     37.7%          91.9%
#
# per_fold and fold_centered collapse at 50% because the large cluster's fold
# holds just two clusters, and centering each fold at its own estimate
# removes most of the large cluster's deviation. Before the leverage (CR3)
# correction, varimpact's column read about 90% at 30% and 70% at 50%.
