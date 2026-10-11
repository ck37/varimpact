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
#' Equals \code{var(ic)} when every cluster is a singleton.
#' @noRd
one_way_var = function(ic, g) {
  G = length(unique(g))
  if (G < 2L) {
    stop("At least 2 clusters are needed to estimate a clustered variance.")
  }
  s = rowsum(ic - mean(ic), g)
  (G / (G - 1)) * sum(s^2) / length(ic)
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
#'
#' @return Numeric scalar.
#' @noRd
cluster_var = function(ic, id = NULL) {
  id = check_cluster_id(id, length(ic))
  if (is.null(id)) return(stats::var(ic))

  k = ncol(id)
  v = 0
  for (m in seq_len(k)) {
    for (cols in utils::combn(k, m, simplify = FALSE)) {
      g = if (length(cols) == 1L) id[[cols]] else
        interaction(id[cols], drop = TRUE)
      v = v + (-1)^(m + 1) * one_way_var(ic, g)
    }
  }

  if (k > 1L && v <= 0) {
    v = max(vapply(id, one_way_var, numeric(1), ic = ic))
    warning("The multiway clustered variance was not positive; using the ",
            "largest one-way clustered variance instead.", call. = FALSE)
  }
  v
}

#' Variance of one fold's influence curve, clustered when ids are available
#'
#' The single entry point used by vim_numerics() and vim_factors(). \code{rows}
#' are the original row numbers of the observations behind \code{ic}, carried
#' through apply_tmle_to_validation() and estimate_pooled_results(), so the
#' cluster ids are looked up by row rather than assumed to line up by position.
#'
#' @param ic Influence curve values for one validation fold.
#' @param rows Row numbers in the full data for each element of \code{ic}.
#' @param cluster_id Output of check_cluster_id() for the full data, or NULL.
#'
#' @return Numeric scalar, or NA when there are too few observations or
#'   clusters in the fold.
#' @noRd
ic_var = function(ic, rows = NULL, cluster_id = NULL) {
  if (length(ic) < 2L) return(NA_real_)
  if (is.null(cluster_id)) return(stats::var(ic))
  if (is.null(rows) || length(rows) != length(ic)) {
    stop("Internal error: the influence curve has ", length(ic),
         " values but ", length(rows), " row indices, so cluster ids cannot ",
         "be matched to it.")
  }
  fold_id = cluster_id[rows, , drop = FALSE]
  too_few = vapply(fold_id, function(x) length(unique(x)) < 2L, logical(1))
  if (any(too_few)) {
    warning("A validation fold has fewer than 2 clusters for ",
            paste(names(fold_id)[too_few], collapse = ", "),
            "; its variance is set to NA.", call. = FALSE)
    return(NA_real_)
  }
  cluster_var(ic, fold_id)
}

#' Degrees of freedom for t-based inference with clustered variances
#'
#' min(G_j) - 1 over clustering variables (Cameron, Gelbach & Miller 2011), or
#' Inf with no clustering, in which case the usual normal critical values are
#' used unchanged.
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
          "'; set folds_by to choose another variable.", call. = FALSE)
  cluster_id[[pick]]
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
