#' Custom printing of the varimpact results.
#'
#' Shows the significant and consistent results by default. If there are no
#' consistent results it shows all results. With cluster-robust standard
#' errors (the \code{id} argument of \code{\link{varimpact}}) it first notes
#' the clustering variables and the largest degrees of freedom.
#'
#' @param x Results object from varimpact.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
print.varimpact = function(x, ...) {
  # Say when the standard errors are clustered, since the tables look the same.
  if (!is.null(x$cluster_id)) {
    clusters = vapply(x$cluster_id, function(g) length(unique(g)), numeric(1))
    cat("Cluster-robust standard errors: ",
        paste0(names(x$cluster_id), " (", clusters, " clusters)",
               collapse = ", "),
        ". Intervals and p-values use a t distribution with the df shown ",
        "for each variable (at most ", x$cluster_df, ").\n\n", sep = "")
  }
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
