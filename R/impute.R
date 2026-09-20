#' Median imputation of a numeric data frame
#'
#' @param data Data frame of numeric columns.
#'
#' @return List with \code{data}, the data frame with each missing value
#'   replaced by its column's median, and \code{medians}, the named vector of
#'   column medians used.
#'
#' @importFrom stats median
#' @noRd
impute_median = function(data) {
  medians = vapply(data, stats::median, numeric(1), na.rm = TRUE)
  for (j in seq_along(medians)) {
    missing = is.na(data[[j]])
    if (any(missing)) {
      data[missing, j] = medians[j]
    }
  }
  list(data = data, medians = medians)
}

#' k-nearest-neighbor imputation of a numeric data frame
#'
#' Reproduces what \code{caret::preProcess(method = "knnImpute")} did: every
#' column is centered and scaled, then each row with a missing value is
#' completed from the \code{k} complete rows nearest to it, by Euclidean
#' distance over the columns it does have, taking the mean of their values in
#' each column it is missing. The output stays centered and scaled.
#'
#' A row missing every column has nothing to match neighbors on. caret stopped
#' on such a row; here it is set to 0, the column mean on the scaled scale,
#' which is what median and zero imputation do in spirit: fill it and move on.
#'
#' @param data Data frame of numeric columns.
#' @param k Number of neighbors, capped at the number of complete rows.
#'
#' @return List with \code{data}, the centered, scaled and imputed data frame,
#'   \code{center} and \code{scale}, the column means and standard deviations
#'   applied, \code{k}, and \code{all_missing}, a logical vector marking the
#'   rows that were missing every column.
#'
#' @importFrom stats sd complete.cases
#' @noRd
impute_knn = function(data, k = 5L) {
  x = as.matrix(data)

  # apply(..., mean) rather than colMeans(): the two sum differently and can
  # disagree in the last bit, and this is what the caret version computed.
  center = apply(x, 2, mean, na.rm = TRUE)
  scale = apply(x, 2, stats::sd, na.rm = TRUE)
  # A constant column, or one with a single observed value, has no usable
  # standard deviation. Leave it centered only; caret did the same.
  scale[is.na(scale) | scale == 0] = 1
  x = sweep(x, 2, center, "-")
  x = sweep(x, 2, scale, "/")

  complete = stats::complete.cases(x)
  all_missing = rowSums(!is.na(x)) == 0

  if (!all(complete)) {
    reference = x[complete, , drop = FALSE]
    if (nrow(reference) == 0) {
      stop("Cannot impute with knn: no row is complete in every numeric covariate.")
    }
    k = min(k, nrow(reference))

    for (i in which(!complete & !all_missing)) {
      observed = which(!is.na(x[i, ]))
      neighbors = RANN::nn2(reference[, observed, drop = FALSE],
                            x[i, observed, drop = FALSE],
                            k = k)$nn.idx
      x[i, -observed] =
        apply(reference[neighbors, -observed, drop = FALSE], 2, mean)
    }
    x[all_missing, ] = 0
  }

  out = data
  out[] = x
  list(data = out, center = center, scale = scale, k = k,
       all_missing = all_missing)
}
