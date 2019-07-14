#' Remove columns from a dataframe if they do not have sufficient variation.
#'
#' @param data Dataframe or matrix
#' @param quantile_probs The probabilities corresponding to the quantiles that will be compared.
#' @param verbose If TRUE output additional details during execution.
#'
#' @return New dataframe with the restriction applied.
#'
#' @seealso quantiles_equivalent
#'
#' @export
restrict_by_quantiles =
  function(data,
           quantile_probs = c(0.1, 0.9),
           verbose = FALSE) {

  # Drop column if the two quantiles have the same sample value (i.e. difference = 0).
  # True = remove, False = keep
  # TODO: support parallelization, e.g. for very wide datasets.
  drop_cols = sapply(1:ncol(data),
                     function(i) quantiles_equivalent(data[, i], quantile_probs))

  if (verbose) {
    cat("restrict_by_quantiles(): dropping", paste(drop_cols, collapse = ", "), "\n")
  }

  # Restrict to variables with sufficient variation.
  data = data[, !drop_cols, drop = FALSE]

  return(data)
}
