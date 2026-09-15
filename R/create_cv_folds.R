#' Stratified CV to insure balance (by one grouping variable, Y)
#'
#' @param V number of folds
#' @param Y Outcome variable. If binary will be used for stratification.
#' @param verbose If T will display extra output.
#'
#' @return Vector of fold assignments.
#'
#' @importFrom cvTools cvFolds
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
      # Record how many observations are in this stratum.
      n = length(rows)
      folds = cvTools::cvFolds(n, K = V, R = 1, type = "random")$which
      out[rows] = folds
    }
    if (verbose) {
      cat("Cross-validation fold breakdown:\n")
      print(table(Y, "Fold"=out, useNA="ifany"))
    }
  } else {
    # More than 2 Ys, so don't stratify.
    xx = cvTools::cvFolds(nn, K = V, R = 1, type = "random")$which
    out = xx
  }
  return(out)
}
