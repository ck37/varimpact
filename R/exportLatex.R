
#' Export varImpact results as Latex tables
#'
#' Outputs results from varImpact() into three Latex tables: consistent results,
#' all results, and per-fold results.
#'
#' Creates three Latex table files:
#'  \itemize{
#'    \item varimpConsistent.tex - the ``consistent'' significant results, meaning those with consistent categories chosen as comparison groups among factors and consistent ordering for numeric variables.
#'    \item varimpAll.tex - the file with cross-validated average variable impacts ordered by statistical significance.
#'    \item varimpByV.tex - the comparison levels used within each validation sample.  Either integer ordering of factors or short-hand for percentile cut-off (0-1 is the 10th percentile, 10+ is the 100th percentile)
#' }
#'
#' @param impact_results Result object from previous varImpact() call.
#' @param outname (Optional) String that is prepended to filenames.
#' @param dir (Optional) Directory to save the results, defaults to current directory.
#' @param digits Digits to round numbers, passed through to xtable.
#'
#' @seealso
#' \code{\link[varImpact]{varImpact}}
#'
#' @export
exportLatex = function(impact_results, outname = "", dir = ".", digits = 4) {
  print(xtable::xtable(impact_results$results_by_fold,
                       caption = "Variable Importance Results By Estimation Sample",
                       label = "byV", digits = digits),
        type = "latex", file = paste0(paste(dir, outname, sep="/"), "varimpByV.tex"),
        caption.placement = "top", include.rownames = T)

  print(xtable::xtable(impact_results$results_all,
                       caption = "Variable Importance Results for Combined Estimates",
                       label = "allRes", digits = digits),
        type = "latex", file = paste0(paste(dir, outname, sep="/"), "varimpAll.tex"),
        caption.placement = "top", include.rownames = T)


  print(xtable::xtable(impact_results$results_consistent,
                       caption = "Subset of of Significant and ``Consistent'' Results",
                       label = "consisRes", digits = digits),
        type = "latex", file = paste0(paste(dir, outname, sep="/"), "varimpConsistent.tex"),
        caption.placement = "top", include.rownames = T)
}
