#' Stratified CV to insure balance (by one grouping variable, Y)
#'
#' @param V number of folds
#' @param Y Outcome variable. If binary will be used for stratification.
#' @param verbose If T will display extra output.
#' @param id Optional grouping vector, one value per observation, from
#'   cluster_fold_groups(). When given, every group is kept within a single
#'   fold. NULL (the default) leaves fold assignment exactly as before.
#'
#' @return Vector of fold assignments.
create_cv_folds = function(V, Y, verbose = F, id = NULL) {
  if (!is.null(id)) {
    return(assign_group_folds(id, V, Y, verbose = verbose))
  }
  # Ignore missing outcomes when deciding whether Y is binary, so that a binary
  # outcome with some missingness is still stratified.
  Ys = unique(Y[!is.na(Y)])
  nys = length(Ys)
  nn = length(Y)
  # Binary outcome so we can do stratified fold generation.
  if (nys == 2) {
    out = rep(NA, nn)
    # Observations with a missing outcome form their own stratum, so that they
    # are spread evenly across the folds rather than clustered in one.
    strata = lapply(Ys, function(y) which(!is.na(Y) & Y == y))
    if (anyNA(Y)) {
      strata = c(strata, list(which(is.na(Y))))
    }
    for (rows in strata) {
      out[rows] = assign_folds(length(rows), V)
    }
    if (verbose) {
      cat("Cross-validation fold breakdown:\n")
      print(table(Y, "Fold"=out, useNA="ifany"))
    }
  } else {
    # More than 2 Ys, so don't stratify.
    out = assign_folds(nn, V)
  }
  return(out)
}

#' Assign observations to folds of as equal size as possible
#'
#' Observation i goes to fold ((i - 1) mod V) + 1, so within a stratum the
#' folds interleave by row order: 1, 2, ..., V, 1, 2, ... This is exactly what
#' \code{cvTools::cvFolds(n, K = V, type = "random")$which} returned, which
#' create_cv_folds() used to read: \code{$which} is the fold of the
#' \emph{permuted} observation, \code{rep(seq_len(K), length.out = n)}, and
#' the permutation itself is in \code{$subsets}, which was never read. So the
#' assignment has always been deterministic, whatever the seed.
#'
#' It is kept that way here on purpose. Drawing the folds at random
#' (\code{rep_len(seq_len(V), n)[sample.int(n)]}) changes the training sets
#' every adjustment step sees, and on the mlbench BreastCancer data that
#' exposed a case where \code{hopach::hopach(mss = "mean")} in
#' reduce_dimensions() never returns. Randomizing the folds therefore needs a
#' guard around HOPACH first, and is a change in its own right rather than
#' part of dropping the cvTools dependency.
#'
#' The unused permutation still consumed n random numbers, and everything
#' downstream that draws from the RNG (SuperLearner's own cross-validation
#' folds, the per-variable future seeds) continues from where it left off. So
#' the same draw is made and discarded here, and a fit from a given seed is
#' identical to one from the previous version. Drop it when the folds become
#' random, since results change then anyway.
#'
#' @param n Number of observations.
#' @param V Number of folds.
#'
#' @return Integer vector of length \code{n} with the fold of each observation.
#'   Fold sizes differ by at most one; when \code{n < V} some folds get no
#'   observation.
#'
#' @noRd
assign_folds = function(n, V) {
  # Consumed and discarded, to keep the RNG stream where cvFolds() left it.
  sample.int(n)
  rep_len(seq_len(V), n)
}

#' Assign whole groups (clusters) to folds
#'
#' Groups are shuffled, then placed largest first into the fold that currently
#' has the fewest observations, which keeps fold sizes close. For a binary
#' outcome, groups containing the rarer value (the events) are placed first,
#' each into the fold with the fewest events so far, so events are spread
#' across folds as evenly as whole groups allow. Missing outcomes count as
#' non-events here. Before either rule, every fold is given two groups when
#' there are enough: otherwise a group holding about 1 / V of the data or
#' more fills a fold by itself, and then a variable measured at the group
#' level has a single value in the other folds' training data and cannot be
#' estimated, and SuperLearner cannot keep groups together within that
#' training data.
#'
#' @param id Grouping vector, one value per observation.
#' @param V Number of folds.
#' @param Y Outcome vector.
#' @param verbose If TRUE, print the fold breakdown.
#'
#' @return Integer vector of length \code{length(Y)} with each observation's
#'   fold.
#' @noRd
assign_group_folds = function(id, V, Y, verbose = FALSE) {
  if (length(id) != length(Y)) {
    stop("The cluster grouping has ", length(id), " values but Y has ",
         length(Y), ".")
  }
  groups = unique(id)
  G = length(groups)
  if (G < V) {
    stop("There are ", G, " clusters to form folds on, fewer than V = ", V,
         ". Reduce V or cluster on a finer variable.")
  }

  group_index = match(id, groups)
  size = tabulate(group_index, nbins = G)
  Ys = unique(Y[!is.na(Y)])
  binary = length(Ys) == 2L
  events = if (binary) {
    # Balance whichever value is rarer: usually Y = 1, but for a missingness
    # indicator, say, it is the 0s that are scarce.
    counts = table(Y[!is.na(Y)])
    rare = as.character(Y) == names(counts)[which.min(counts)]
    tabulate(group_index[!is.na(Y) & rare], nbins = G)
  } else {
    rep(0L, G)
  }

  # Shuffle so that ties are broken at random, then order: groups with more
  # events first, then larger groups first.
  shuffled = sample.int(G)
  ord = shuffled[order(-events[shuffled], -size[shuffled])]

  fold_n = rep(0L, V)
  fold_events = rep(0L, V)
  fold_groups = rep(0L, V)
  group_fold = integer(G)
  for (g in ord) {
    candidates = seq_len(V)
    # Until every fold has two groups, only folds with the fewest groups are
    # eligible.
    if (min(fold_groups) < 2L) {
      candidates = which(fold_groups == min(fold_groups))
    }
    if (events[g] > 0L) {
      fewest = min(fold_events[candidates])
      candidates = candidates[fold_events[candidates] == fewest]
    }
    f = candidates[which.min(fold_n[candidates])]
    group_fold[g] = f
    fold_n[f] = fold_n[f] + size[g]
    fold_events[f] = fold_events[f] + events[g]
    fold_groups[f] = fold_groups[f] + 1L
  }

  out = group_fold[group_index]
  if (verbose) {
    cat("Cluster-grouped cross-validation:", G, "clusters in", V, "folds.\n")
    if (binary) {
      print(table(Y, "Fold" = out, useNA = "ifany"))
    } else {
      print(table("Fold" = out))
    }
  }
  out
}
