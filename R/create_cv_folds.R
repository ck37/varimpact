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
  Ys = unique(Y)
  nys = length(Ys)
  nn = length(Y)
  # Binary outcome so we can do stratified fold generation.
  if (nys == 2) {
    out = rep(NA, nn)
    for (i in 1:nys) {
      # Record how many observations have this Y value.
      n = sum(Y == Ys[i])
      folds = cvTools::cvFolds(n, K = V, R = 1, type = "random")$which
      #if (verbose) {
      #  cat("Y", i, "is", Ys[i], "count:", sum(Y == Ys[i]), "n=", n, "fold length:",
      #      length(folds), "\n")
      #}
      out[Y == Ys[i]] = folds
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
