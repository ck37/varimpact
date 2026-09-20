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

#' Randomly assign observations to folds of as equal size as possible
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
  # Indexing by sample.int() rather than calling sample() on the vector: for
  # n = 1 the vector has one element and sample() would treat it as a range.
  rep_len(seq_len(V), n)[sample.int(n)]
}
