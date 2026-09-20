#' Stratified CV to insure balance (by one grouping variable, Y)
#'
#' @param V number of folds
#' @param Y Outcome variable. If binary will be used for stratification.
#' @param verbose If T will display extra output.
#'
#' @return Vector of fold assignments.
create_cv_folds = function(V, Y, verbose = F) {
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
