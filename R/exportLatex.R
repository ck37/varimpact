
#' Export varimpact results as Latex tables
#'
#' Outputs results from varimpact() into three Latex tables: consistent results,
#' all results, and per-fold results.
#'
#' Creates three Latex table files:
#'  \itemize{
#'    \item varimpConsistent.tex - the ``consistent'' significant results, meaning those with consistent categories chosen as comparison groups among factors and consistent ordering for numeric variables.
#'    \item varimpAll.tex - the file with cross-validated average variable impacts ordered by statistical significance.
#'    \item varimpByV.tex - the comparison levels used within each validation sample.  Either integer ordering of factors or short-hand for percentile cut-off (0-1 is the 10th percentile, 10+ is the 100th percentile)
#' }
#'
#' @param impact_results Result object from previous varimpact() call.
#' @param outname (Optional) String that is prepended to filenames.
#' @param dir (Optional) Directory to save the results, defaults to current directory.
#' @param digits Digits to round numbers, passed through to xtable.
#' @param ... Additional parameters passed to print.xtable().
#'
#' @seealso
#' \code{\link[varimpact]{varimpact}}
#'
#' @export
# TODO: document return object.
exportLatex = function(impact_results, outname = "", dir = ".", digits = 4, ...) {

  table_byfold = cbind("Variable" = rownames(impact_results$results_by_fold),
                    impact_results$results_by_fold)

  xtable_byfold = xtable::xtable(table_byfold,
            caption = "Variable Importance Results By Estimation Sample",
            label = "byFold",
            digits = digits)

  print(xtable_byfold,
        type = "latex",
        file = paste0(dir, "/", outname, "varimpByFold.tex"),
        caption.placement = "top",
        include.rownames = F,
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

  table_all = cbind("Rank" = 1:nrow(impact_results$results_all),
                    "Variable" = rownames(impact_results$results_all),
                    impact_results$results_all)

  xtable_all = xtable::xtable(table_all,
                   caption = "Variable Importance Results For Combined Estimates",
                   label = "allRes",
                   digits = digits)

  print(xtable_all,
          type = "latex",
          file = paste0(dir, "/", outname, "varimpAll.tex"),
          caption.placement = "top",
          include.rownames = F,
          hline.after = hline.after,
          ...)

  if (nrow(impact_results$results_consistent) > 0) {
    table_consistent = cbind("Rank" = 1:nrow(impact_results$results_consistent),
                           "Variable" = rownames(impact_results$results_consistent),
                           impact_results$results_consistent)
    xtable_consistent = xtable::xtable(table_consistent,
              caption = "Subset of of Significant and ``Consistent'' Results",
              label = "consisRes",
              digits = digits)
  } else {
    # Create a blank dataframe.
    table_consistent = data.frame()
    # Create a blank xtable.
    xtable_consistent = NULL
  }

  print(xtable_consistent,
        type = "latex",
        file = paste0(dir, "/", outname, "varimpConsistent.tex"),
        caption.placement = "top",
        include.rownames = F,
        ...)

  # Return a list with the output results.
  results = list(tables = list(
                   consistent = table_consistent,
                   all = table_all,
                   byfold = table_byfold
                 ),
                 xtables = list(
                   consistent = xtable_consistent,
                   all = xtable_all,
                   byfold = xtable_byfold
                 ))

  return(invisible(results))
}
