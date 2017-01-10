#' Custom printing of the varImpact results.
#'
#' Shows the significant and consistent results by default. If there are no
#' consistent results it shows all results.
#'
#' @param x Results object from varImpact.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
print.varImpact = function(x, ...) {
  # Just print the significant and consistent results.
  if (!is.null(x$results_consistent) && nrow(x$results_consistent) > 0) {
    cat("Significant and consistent results:\n")
    print(x$results_consistent)
  } else if (!is.null(x$results_all)) {
    cat("No significant and consistent results.\n")
    cat("All results:\n")
    print(x$results_all)
  } else {
    cat("No results could be calculated.\n")
  }
  invisible(x)
}
