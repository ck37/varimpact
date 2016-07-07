#' Find the maximum of the squared values.
#'
#' @param x Numeric vector
#'
#' @return Maximum of the squared values, or -Inf if all elements are NA.
#' @export
max_sqr = function(x) {
  # Handle missing data manually so we don't get warnings when it occurs.
  x = na.omit(x)
  if (length(x) == 0) {
    -Inf
  } else {
    max(x^2)
  }
}
