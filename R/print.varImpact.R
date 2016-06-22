#' Custom printing of the varImpact results.
#'
#' Shows the consistent results by default. If there are no consistent results
#' it shows all results.
#'
#' @param x Results object from varImpact.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
print.varImpact = function(x, ...) {
  # Just print the consistent results.
  if (nrow(x$results_consistent) > 0) {
    cat("Consistent results:\n")
    print(x$results_consistent)
  } else {
    cat("No consistent results.\n")
    cat("All results:\n")
    print(x$results_all)
  }
  invisible(x)
}
