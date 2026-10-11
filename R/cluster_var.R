# Cluster-robust variance estimation for influence curves.
#
# varimpact's standard errors come from the variance of the influence curve
# (IC) on each validation fold: var(IC) / n assumes every row is an independent
# draw. When rows are correlated within clusters (siblings in a family, patients
# in a hospital, students in a school), that understates the variance. The
# cluster-robust version sums the IC within each cluster and treats the cluster
# totals as the independent units, which is the influence-curve analogue of the
# sandwich estimator in sandwich::vcovCL().
#
# With several clustering variables the Cameron, Gelbach & Miller (2011)
# inclusion-exclusion estimator is used, e.g. V = V_A + V_B - V_(A x B).
#
# References:
#   Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). Robust inference
#     with multiway clustering. Journal of Business & Economic Statistics,
#     29(2), 238-249.
#   Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to
#     cluster-robust inference. Journal of Human Resources, 50(2), 317-372.

#' Validate and standardize cluster identifiers
#'
#' Accepts NULL, a single vector, or a data frame / matrix with one column per
#' clustering variable, and returns a data frame of integer codes. Two kinds of
#' column are dropped because they cannot change the variance:
#' \itemize{
#'   \item a column in which every cluster has one observation (clustering on
#'     singletons is the same as not clustering), and
#'   \item a column nested inside another column, e.g. family inside site:
#'     clustering on the outer variable already sums the inner clusters
#'     together, and the multiway formula reduces exactly to the outer one-way
#'     variance.
#' }
#'
#' @param id NULL, a vector of length \code{n}, or a data frame or matrix with
#'   \code{n} rows.
#' @param n Number of observations.
#' @param verbose If TRUE, report columns that are dropped.
#'
#' @return NULL when there is nothing to cluster on, otherwise a data frame of
#'   integer-coded cluster ids with at least one column.
#' @noRd
check_cluster_id = function(id, n, verbose = FALSE) {
  if (is.null(id)) return(NULL)
  if (is.null(dim(id))) {
    if (is.list(id)) id = as.data.frame(id, stringsAsFactors = FALSE)
    else id = data.frame(id = id, stringsAsFactors = FALSE)
  }
  id = as.data.frame(id, stringsAsFactors = FALSE)
  if (ncol(id) == 0L) return(NULL)
  if (nrow(id) != n) {
    stop("The cluster id has ", nrow(id), " rows but there are ", n,
         " observations.")
  }
  if (anyNA(id)) {
    stop("The cluster id contains missing values; every observation needs a ",
         "cluster.")
  }
  if (is.null(names(id)) || any(names(id) == "") || anyDuplicated(names(id))) {
    names(id) = paste0("id", seq_len(ncol(id)))
  }
  id[] = lapply(id, function(x) as.integer(factor(x)))

  # Singleton-only columns carry no clustering.
  singleton = vapply(id, function(x) length(unique(x)) == n, logical(1))
  if (verbose && any(singleton)) {
    message("Cluster variable(s) ", paste(names(id)[singleton], collapse = ", "),
            " have one observation per cluster and are ignored.")
  }
  id = id[, !singleton, drop = FALSE]
  if (ncol(id) == 0L) return(NULL)

  # Drop columns nested within another column (this also removes duplicates).
  k = ncol(id)
  keep = rep(TRUE, k)
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i != j && keep[i] && keep[j] && is_nested(id[[i]], id[[j]])) {
        if (verbose) {
          message("Cluster variable '", names(id)[i], "' is nested in '",
                  names(id)[j], "'; clustering on '", names(id)[j],
                  "' already accounts for it, so it is not used separately.")
        }
        keep[i] = FALSE
      }
    }
  }
  id[, keep, drop = FALSE]
}

#' TRUE if every level of \code{inner} falls within a single level of
#' \code{outer}.
#' @noRd
is_nested = function(inner, outer) {
  nrow(unique(data.frame(inner, outer))) == length(unique(inner))
}

#' One-way clustered variance of an influence curve, per observation.
#'
#' Equals \code{var(ic)} when every cluster is a singleton. Without weights
#' this is the CR1 version, with a G / (G - 1) factor. With \code{weights},
#' each row's weight in the estimate, it is the jackknife-style CR3 version:
#' each cluster's sum is divided by 1 - h, where h is the cluster's share of
#' the weights, with a (G - 1) / G factor. The two agree when every cluster
#' has the same share.
#' @noRd
one_way_var = function(ic, g, weights = NULL) {
  G = length(unique(g))
  if (G < 2L) {
    stop("At least 2 clusters are needed to estimate a clustered variance.")
  }
  s = rowsum(ic - mean(ic), g)
  if (is.null(weights)) {
    return((G / (G - 1)) * sum(s^2) / length(ic))
  }
  # Jackknife-style (CR3) version: each cluster's sum is inflated by its
  # leverage, its share of the weights, so a dominant cluster is not
  # understated. With equal shares (1 / G) this equals the G / (G - 1) version.
  h = rowsum(weights, g) / sum(weights)
  ((G - 1) / G) * sum((s / (1 - h))^2) / length(ic)
}

#' Cluster-robust variance of an influence curve
#'
#' Returns a per-observation variance, so \code{cluster_var(ic, id) / n}
#' estimates the variance of the estimator, matching how \code{var(ic)} is used
#' elsewhere in varimpact. With no id, or with every observation in its own
#' cluster, the result equals \code{var(ic)} exactly.
#'
#' With more than one clustering variable the Cameron, Gelbach & Miller (2011)
#' inclusion-exclusion estimator is used: for two variables,
#' \eqn{V = V_A + V_B - V_{A \cap B}}, where \eqn{V_{A \cap B}} clusters on each
#' unique combination of A and B. Each term carries its own \eqn{G/(G-1)}
#' small-sample factor, as in \code{sandwich::vcovCL(multi0 = FALSE)}. If the
#' combination is not positive, which can happen in finite samples, the
#' largest one-way variance is returned with a warning.
#'
#' @param ic Numeric influence curve values.
#' @param id NULL, a vector, or a data frame of cluster ids aligned with
#'   \code{ic}.
#' @param weights NULL for the G / (G - 1) (CR1) version, or each row's weight
#'   in the estimate for the jackknife-style CR3 version (see one_way_var()).
#'
#' @return Numeric scalar.
#' @noRd
cluster_var = function(ic, id = NULL, weights = NULL) {
  id = check_cluster_id(id, length(ic))
  if (is.null(id)) return(stats::var(ic))

  k = ncol(id)
  v = 0
  for (m in seq_len(k)) {
    for (cols in utils::combn(k, m, simplify = FALSE)) {
      g = if (length(cols) == 1L) id[[cols]] else
        interaction(id[cols], drop = TRUE)
      v = v + (-1)^(m + 1) * one_way_var(ic, g, weights)
    }
  }

  if (k > 1L && v <= 0) {
    v = max(vapply(id, one_way_var, numeric(1), ic = ic, weights = weights))
    warning("The multiway clustered variance was not positive; using the ",
            "largest one-way clustered variance instead.", call. = FALSE)
  }
  v
}

#' Variance of one fold's influence curve, clustered when ids are available
#'
#' Used for the per-fold variances that vim_numerics() and vim_factors() report
#' alongside each fold's estimate. With clustering, the standard errors in the
#' final results come from cv_cluster_var() instead. \code{rows} are the
#' original row numbers of the observations behind \code{ic}, carried through
#' apply_tmle_to_validation() and estimate_pooled_results(), so the cluster ids
#' are looked up by row rather than assumed to line up by position.
#'
#' @param ic Influence curve values for one validation fold.
#' @param rows Row numbers in the full data for each element of \code{ic}.
#' @param cluster_id Output of check_cluster_id() for the full data, or NULL.
#'
#' @return Numeric scalar, or NA when there are too few observations or
#'   clusters in the fold. A fold can hold a single cluster (one large hospital,
#'   say), and a variance within it cannot be estimated, but that does not
#'   affect the clustered standard error, so this is not a warning.
#' @noRd
ic_var = function(ic, rows = NULL, cluster_id = NULL) {
  if (length(ic) < 2L) return(NA_real_)
  if (is.null(cluster_id)) return(stats::var(ic))
  check_ic_rows(ic, rows)
  fold_id = cluster_id[rows, , drop = FALSE]
  if (any(vapply(fold_id, function(x) length(unique(x)) < 2L, logical(1)))) {
    return(NA_real_)
  }
  cluster_var(ic, fold_id)
}

#' Cluster-robust variance of the cross-validated estimate
#'
#' varimpact's estimate is the average of the V fold estimates, so to first
#' order its error is \eqn{\sum_v \frac{1}{V n_v} \sum_{i \in v} D_i}, where
#' \eqn{D_i} is observation i's influence curve (Zheng & van der Laan 2010,
#' Asymptotic theory for cross-validated targeted maximum likelihood
#' estimation, U.C. Berkeley Division of Biostatistics Working Paper 273; and,
#' for parameters whose levels are chosen within each training fold as
#' varimpact's are, Hubbard, Kherad-Pajouh & van der Laan 2016, Statistical
#' inference for data adaptive target parameters, International Journal of
#' Biostatistics 12(1)). With clusters as the independent units and each
#' cluster within one fold, that sum splits into independent cluster totals,
#' and its variance is the sum of their squares. The tmle package takes the
#' same influence-curve approach with its id argument, averaging rather than
#' summing within id; the two agree for clusters of equal size, and summing
#' matches varimpact's estimate, which weights observations equally. This
#' scales each fold's curve
#' by \eqn{n / (V n_v)}, stacks the folds, and clusters once over all n
#' observations. Unlike averaging per-fold clustered variances, this counts a
#' cluster that spans several folds (as crossed clustering variables can) as
#' one cluster, and it does not need two or more clusters in every fold.
#'
#' The curves must be centered at the overall estimate, not at each fold's own
#' estimate: see cv_cluster_vars(). Centering within each fold would remove the
#' between-cluster variation of any fold that holds a single cluster, such as
#' one large hospital, and understate the variance badly.
#'
#' Each cluster's sum is scaled by 1 / (1 - h), where h is the cluster's share
#' of the weight in the average, with a (G - 1) / G factor: the jackknife-style
#' CR3 variance (MacKinnon, Nielsen & Webb 2023, Cluster-robust inference: A
#' guide to empirical practice, Journal of Econometrics 232(2)), as in
#' sandwich::vcovCL(type = "HC3", cadjust = FALSE). For the average of fold
#' means it equals the cluster jackknife, and with clusters of equal weight it
#' equals the usual G / (G - 1) (CR1) version. It corrects the
#' understatement from a cluster that holds a large share of the data: in
#' simulations of varimpact() with one cluster holding half the rows, false
#' positives for a null variable fell from 7.6% to 4.6%.
#'
#' @param ic_list List with one influence-curve vector per fold.
#' @param rows_list List with the matching row numbers in the full data.
#' @param cluster_id Output of check_cluster_id() for the full data.
#'
#' @return A per-observation variance, so that dividing by n gives the variance
#'   of the estimate, like the average of the per-fold variances it replaces.
#'   NA when any fold has fewer than 2 observations, matching the unclustered
#'   path, where such a fold's variance is NA.
#' @noRd
cv_cluster_var = function(ic_list, rows_list, cluster_id) {
  sizes = lengths(ic_list)
  if (length(sizes) == 0L || any(sizes < 2L)) return(NA_real_)
  for (v in seq_along(ic_list)) check_ic_rows(ic_list[[v]], rows_list[[v]])
  n = sum(sizes)
  V = length(ic_list)
  scaled = unlist(lapply(seq_len(V), function(v) {
    ic_list[[v]] * n / (V * sizes[v])
  }))
  # Each row's weight in the average of the fold estimates, for the leverage
  # (CR3) correction.
  weights = rep(1 / (V * sizes), sizes)
  rows = unlist(rows_list)
  cluster_var(scaled, cluster_id[rows, , drop = FALSE], weights)
}

#' Clustered variances of the risk difference and the log relative risk
#'
#' estimate_pooled_results() centers each fold's influence curve at that
#' fold's estimate. Shifting it by the fold's deviation from the average over
#' folds re-centers it at the estimate varimpact reports, so the clustered sums
#' in cv_cluster_var() keep the variation between folds. The two contrasts are
#' then formed as vim_numerics() and vim_factors() form their per-fold
#' variances. Both levels are predicted on the same validation rows, so the
#' high level's row numbers serve for both.
#'
#' @param pooled_min,pooled_max Output of estimate_pooled_results() for the
#'   low and high levels.
#' @param V Number of folds.
#' @param cluster_id Output of check_cluster_id() for the full data.
#'
#' @return List with elements \code{rd} and \code{log_rr}.
#' @noRd
cv_cluster_vars = function(pooled_min, pooled_max, V, cluster_id) {
  recenter = function(pooled) {
    lapply(seq_len(V), function(v) {
      pooled$influence_curves[[v]] + pooled$thetas[[v]] - mean(pooled$thetas)
    })
  }
  ic_min = recenter(pooled_min)
  ic_max = recenter(pooled_max)
  rd = lapply(seq_len(V), function(v) ic_max[[v]] - ic_min[[v]])
  log_rr = lapply(seq_len(V), function(v) {
    ic_max[[v]] / pooled_max$thetas[[v]] - ic_min[[v]] / pooled_min$thetas[[v]]
  })
  list(rd = cv_cluster_var(rd, pooled_max$rows, cluster_id),
       log_rr = cv_cluster_var(log_rr, pooled_max$rows, cluster_id))
}

#' Degrees of freedom for one variable's contrast
#'
#' The global degrees of freedom count every cluster, but a variable's estimate
#' compares its low and high level, and when only a few clusters hold one of
#' those levels (a hospital-level policy in a handful of hospitals, say) the
#' clustered variance of the contrast rests on those few clusters. Normal or
#' t(G - 1) critical values then over-reject badly (Imbens & Kolesar 2016,
#' Robust standard errors in small samples, Review of Economics and Statistics
#' 98(4); MacKinnon & Webb 2017). As a simple, conservative stand-in for the
#' Bell-McCaffrey effective degrees of freedom (Bell & McCaffrey 2002, Bias
#' reduction in standard errors for linear regression with multi-stage
#' samples, Survey Methodology 28(2)), which is not itself a published rule,
#' this counts, in each fold's
#' observed validation rows, the clusters holding the low level and those
#' holding the high level, takes the fewer, sums over folds, and subtracts
#' one. Counting fold by fold matters because folds can choose different
#' levels as low and high. compile_results() caps the result at the global df.
#' With several clustering variables the smallest count is used. In
#' simulations this kept false positives for cluster-level factors nearer 5%
#' than a single global df or a Welch-type df that also counts the other
#' level's clusters, at the cost of being conservative with very few clusters.
#'
#' @param fold_results List with one element per fold, each holding
#'   \code{level_min$val_preds} and \code{level_max$val_preds} from
#'   apply_tmle_to_validation().
#' @param cluster_id Output of check_cluster_id() for the full data.
#'
#' @return Numeric scalar.
#' @noRd
contrast_df = function(fold_results, cluster_id) {
  # Clusters holding a level among one fold's observed validation rows, for
  # each clustering variable.
  clusters = function(preds) {
    rows = preds$row_index[preds$A == 1 & preds$delta == 1]
    vapply(cluster_id, function(g) length(unique(g[rows])), numeric(1))
  }
  per_fold = lapply(fold_results, function(x) {
    pmin(clusters(x$level_min$val_preds), clusters(x$level_max$val_preds))
  })
  min(Reduce(`+`, per_fold)) - 1
}

#' Stop if an influence curve and its row numbers do not line up.
#' @noRd
check_ic_rows = function(ic, rows) {
  if (is.null(rows) || length(rows) != length(ic)) {
    stop("Internal error: the influence curve has ", length(ic),
         " values but ", length(rows), " row indices, so cluster ids cannot ",
         "be matched to it.")
  }
}

#' Degrees of freedom for t-based inference with clustered variances
#'
#' min(G_j) - 1 over clustering variables (Cameron, Gelbach & Miller 2011), or
#' Inf with no clustering, in which case the usual normal critical values are
#' used unchanged. t(G - 1) critical values are the usual reference with
#' clustered standard errors (Cameron & Miller 2015; MacKinnon, Nielsen & Webb
#' 2023).
#'
#' @param cluster_id Output of check_cluster_id(), or NULL.
#' @noRd
cluster_df = function(cluster_id) {
  if (is.null(cluster_id)) return(Inf)
  min(vapply(cluster_id, function(x) length(unique(x)), numeric(1))) - 1
}

#' Collapse the cluster ids into one grouping for CV fold assignment
#'
#' CV-TMLE needs each validation fold to be independent of its training
#' folds, so every cluster has to fall entirely inside one fold. With several
#' clustering variables, the grouping that keeps all of them intact is the set
#' of connected components: two rows are linked when they share a cluster on
#' any variable. For nested variables this is the outer variable. Crossed
#' variables can chain nearly everything into one component; when fewer than
#' V components remain, the coarsest variable that still has at least V
#' clusters is used instead, with a warning, or the variable named in
#' \code{folds_by}.
#'
#' @param cluster_id Output of check_cluster_id(), or NULL.
#' @param V Number of CV folds.
#' @param folds_by Optional name of the cluster variable to form folds on.
#'
#' @return NULL, or an integer vector with one grouping value per row.
#' @noRd
cluster_fold_groups = function(cluster_id, V, folds_by = NULL) {
  if (is.null(cluster_id)) {
    if (!is.null(folds_by)) {
      stop("folds_by was given but there is no cluster id to form folds on.")
    }
    return(NULL)
  }
  if (!is.null(folds_by)) {
    if (!folds_by %in% names(cluster_id)) {
      stop("folds_by = '", folds_by, "' is not one of the cluster variables ",
           "in use: ", paste(names(cluster_id), collapse = ", "), ". Nested ",
           "and singleton cluster variables are dropped before folding.")
    }
    return(cluster_id[[folds_by]])
  }
  if (ncol(cluster_id) == 1L) return(cluster_id[[1]])

  # Connected components by label propagation across the clustering variables.
  comp = seq_len(nrow(cluster_id))
  repeat {
    old = comp
    for (col in cluster_id) comp = stats::ave(comp, col, FUN = min)
    if (identical(comp, old)) break
  }
  comp = as.integer(factor(comp))
  if (length(unique(comp)) >= V) return(comp)

  G = vapply(cluster_id, function(x) length(unique(x)), numeric(1))
  ok = which(G >= V)
  if (length(ok) == 0L) {
    stop("No cluster variable has at least V = ", V, " clusters; reduce V.")
  }
  pick = ok[which.min(G[ok])]
  warning("The cluster variables are crossed, so they cannot all be kept ",
          "within single folds. Folding on '", names(cluster_id)[pick],
          "': training and validation folds then share clusters of ",
          paste0("'", names(cluster_id)[-pick], "'", collapse = ", "),
          ", which can bias estimates for variables that vary with them. Set ",
          "folds_by to the variable that the variables of interest vary with.",
          call. = FALSE)
  cluster_id[[pick]]
}

#' Drop predictors that cannot be compared across cluster-grouped folds
#'
#' With whole clusters assigned to folds, a factor whose every level falls
#' within a single fold group (the cluster id itself, or a finer label such as
#' a person within a family) never has the same level in a training fold and
#' its validation fold, so its importance cannot be estimated. Left in, such a
#' variable dropped out of the results with only "Failed in
#' results_by_level()" printed. Numeric columns stay even when they are
#' constant within clusters (hospital size, say), since they are binned into
#' ranges that span clusters.
#'
#' @param data Data frame of predictors.
#' @param fold_groups Output of cluster_fold_groups(), or NULL.
#'
#' @return \code{data} without those columns, with a message naming them.
#' @noRd
drop_cluster_level_columns = function(data, fold_groups) {
  if (is.null(fold_groups)) return(data)
  # Column by column rather than vapply(data, ...), so a matrix works too.
  drop = vapply(seq_len(ncol(data)), function(j) {
    x = data[, j]
    (is.factor(x) || is.character(x)) && is_nested(x, fold_groups)
  }, logical(1))
  if (any(drop)) {
    message("Not analyzing ", paste(colnames(data)[drop], collapse = ", "),
            ": each value falls within a single cluster, so with whole ",
            "clusters assigned to CV folds it never appears in both a ",
            "training fold and its validation fold.")
    data = data[, !drop, drop = FALSE]
    if (ncol(data) == 0L) {
      stop("data has no variables left that vary across clusters.")
    }
  }
  data
}

#' Number of SuperLearner CV folds to use with a (possibly clustered) id
#'
#' SuperLearner::CVFolds() assigns whole ids to folds and fails when there are
#' fewer distinct ids than folds. Returns \code{V} unchanged when there are
#' enough ids, otherwise the number of distinct ids (possibly < 2, which the
#' caller must handle by dropping the id).
#' @noRd
cluster_sl_folds = function(id, V) {
  as.integer(min(V, length(unique(id))))
}
