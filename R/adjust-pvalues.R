#' Holm and Benjamini-Hochberg adjusted p-values
#'
#' Wraps \code{stats::p.adjust()} to reproduce what
#' \code{multtest::mt.rawp2adjp(p, c("Holm", "BH"))} returned: a matrix with
#' the raw p-values and the two adjusted versions, one row per input.
#'
#' @param p Vector of raw p-values. Missing values are allowed.
#'
#' @return Numeric matrix with columns \code{rawp}, \code{Holm} and
#'   \code{BH}, in the order of \code{p}. Rows whose raw p-value is missing
#'   have missing adjusted values.
#'
#' @importFrom stats p.adjust
#' @noRd
adjust_pvalues = function(p) {
  # multtest counted every p-value toward the number of tests, missing ones
  # included (its na.rm = FALSE default). p.adjust() would count only the
  # non-missing ones, so pass the full count to keep the adjusted values the
  # same as before. A missing p-value belongs to a variable whose estimate
  # failed, and compile_results() drops it from the output afterward.
  n = length(p)
  adjusted = cbind(rawp = p,
                   Holm = stats::p.adjust(p, method = "holm", n = n),
                   BH = stats::p.adjust(p, method = "BH", n = n))
  adjusted[is.na(p), c("Holm", "BH")] = NA_real_
  adjusted
}
