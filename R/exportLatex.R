
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
#' @param ... Additional parameters passed to print.xtable().
#'
#' @seealso
#' \code{\link[varImpact]{varImpact}}
#'
#' @export
exportLatex = function(impact_results, outname = "", dir = ".", digits = 4, ...) {
  print(xtable::xtable(impact_results$results_by_fold,
            caption = "Variable Importance Results By Estimation Sample",
            label = "byV",
            digits = digits),
        type = "latex",
        file = paste0(paste(dir, outname, sep="/"), "varimpByFold.tex"),
        #caption.placement = "top",
        include.rownames = T,
        ...)


  # Use hline.after to add a line after the p = 0.05 cut-off.
  signif_cutoff = which(impact_results$results_all[, "Adj. p-value"] > 0.05)
  if (length(signif_cutoff) > 0) {
    signif_row = min(signif_cutoff) - 1
    hline.after = c(-1,0, signif_row, nrow(impact_results$results_all))
  } else {
    # All variables are important.
    hline.after = NULL
  }

  table_all = cbind("Rank"=1:nrow(impact_results$results_all),
                    "Variable"=rownames(impact_results$results_all),
                    impact_results$results_all)

  print(xtable::xtable(table_all,
                  caption = "Variable Importance Results for Combined Estimates",
                  label = "allRes",
                  digits = digits),
          type = "latex",
          file = paste0(paste(dir, outname, sep="/"), "varimpAll.tex"),
          caption.placement = "top",
          include.rownames = F,
          hline.after = hline.after,
          ...)

  if (nrow(impact_results$results_consistent) > 0) {
    consistent_table = cbind("Rank"=1:nrow(impact_results$results_consistent),
                           "Variable"=rownames(impact_results$results_consistent),
                           impact_results$results_consistent)
    consistent_xtable = xtable::xtable(consistent_table,
              caption = "Subset of of Significant and ``Consistent'' Results",
              label = "consisRes",
              digits = digits)
  } else {
    # Create a blank dataframe.
    consistent_table = data.frame()
    # Create a blank xtable.
    consistent_xtable = NULL
  }

  print(consistent_xtable,
        type = "latex",
        file = paste0(paste(dir, outname, sep="/"), "varimpConsistent.tex"),
        caption.placement = "top",
        include.rownames = F,
        ...)

  # Give a default return value, silently.
  # TODO: return a list with the output results.
  invisible(T)
}
