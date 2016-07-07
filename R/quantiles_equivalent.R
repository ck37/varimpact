#' Checks if two quantiles have the same sample value for a given vector.
#' This indicates that the vector has little variation.
#'
#' @param x Data vector. If a factor it is converted to numeric using unclass().
#' @param quantile_probs Vector with two probabilities that specify the quantiles.
#'
#' @return True if the two quantiles are equal, indicating a lack of variation in the sample data.
#'
#' @seealso restrict_by_quantiles
#'
#' @export
quantiles_equivalent = function(x, quantile_probs = c(0.1, 0.9)) {
  if (length(quantile_probs) != 2) {
    warning("Quantiles_equivalent() expects quantile_probs to be a 2-element vector.")
  }
  if (class(x) == "factor") {
    x = unclass(x)
  }
  quantiles = quantile(x, probs = quantile_probs, na.rm = T)
  # Returns True if there is no differences between the first and second quantiles.
  (quantiles[2] - quantiles[1]) == 0
}
