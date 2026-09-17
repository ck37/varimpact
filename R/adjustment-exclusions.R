# Helpers for the adjustment_exclusions argument of varimpact().

#' Adjustment matrix columns belonging to one input variable
#'
#' varimpact() expands the analysis data before it becomes an adjustment
#' matrix, so a single input variable can appear in W under several different
#' names: a numeric \code{V1} stays \code{"V1"}, a factor \code{F1} becomes one
#' indicator column per level (\code{"F1XXb"}, \code{"F1XXc"} - see
#' \code{\link{factors_to_indicators}}), and either one gains \code{"Imiss_V1"}
#' when it has missing values. Matching an exclusion request against
#' \code{colnames(W)} literally would therefore drop a numeric variable but
#' silently leave every indicator of a factor in place.
#'
#' @param name Name of one column of the analysis data.
#' @param w_names Column names of the adjustment matrix.
#'
#' @return Character vector of the columns in \code{w_names} that came from
#'   \code{name}; empty if it contributed none.
#'
#' @seealso \code{\link{exclude_adjustment_vars}}
adjustment_columns = function(name, w_names) {
  w_names[w_names == name |
          startsWith(w_names, paste0(name, "XX")) |
          w_names == paste0("Imiss_", name)]
}

#' Remove excluded adjustment variables from an adjustment matrix
#'
#' Drops the adjustment variables the user excluded for this candidate
#' variable, resolving each requested name to every column it contributed via
#' \code{\link{adjustment_columns}}.
#'
#' This runs inside the per-variable, per-fold loop, so it stays silent: the
#' spelling of \code{adjustment_exclusions} is checked once up front by
#' \code{\link{check_adjustment_exclusions}} rather than warning on every fold.
#'
#' @param W Adjustment matrix (data frame) for the current candidate variable.
#' @param name Name of the candidate variable whose importance is being
#'   estimated.
#' @param adjustment_exclusions Named list as passed to \code{\link{varimpact}}.
#' @param verbose If TRUE, report which columns were dropped.
#'
#' @return \code{W} without the excluded columns, unchanged when nothing was
#'   requested for \code{name} or nothing matched.
#'
#' @seealso \code{\link{varimpact}}
exclude_adjustment_vars = function(W, name, adjustment_exclusions,
                                   verbose = FALSE) {
  requested = adjustment_exclusions[[name]]

  if (is.null(requested) || length(requested) == 0L) {
    return(W)
  }

  drop = unique(unlist(lapply(requested, adjustment_columns,
                              w_names = colnames(W))))

  if (length(drop) == 0L) {
    return(W)
  }

  if (verbose) {
    cat("Excluding", length(drop), "adjustment column(s) for", name, ":",
        paste(drop, collapse = ", "), "\n")
  }

  W[, setdiff(colnames(W), drop), drop = FALSE]
}

#' Validate an adjustment_exclusions list
#'
#' Checks \code{adjustment_exclusions} once against the columns of the analysis
#' data, before any folds are run.
#'
#' A name that matches nothing is almost always a typo, and excluding nothing
#' silently is the dangerous failure here: the user would believe they had
#' adjusted away a mediator or a collider when they had not. Unknown names warn
#' rather than stop, so that a stale entry left over from an earlier version of
#' the data does not abort a long run; only a malformed argument is an error.
#'
#' @param adjustment_exclusions Named list as passed to \code{\link{varimpact}}.
#'   An empty list is accepted and checks nothing.
#' @param data_names Column names of the analysis data.
#'
#' @return Invisibly NULL; called for its warnings and errors.
#'
#' @seealso \code{\link{exclude_adjustment_vars}}
check_adjustment_exclusions = function(adjustment_exclusions, data_names) {
  if (length(adjustment_exclusions) == 0L) {
    return(invisible(NULL))
  }

  if (!is.list(adjustment_exclusions) || is.null(names(adjustment_exclusions))) {
    stop("adjustment_exclusions must be a named list, e.g. ",
         "list(V1 = c(\"V2\", \"V3\")). Each name is a variable whose ",
         "importance is being estimated; each element lists the adjustment ",
         "variables to exclude for it.")
  }

  unknown_targets = setdiff(names(adjustment_exclusions), data_names)
  if (length(unknown_targets) > 0L) {
    warning("adjustment_exclusions names variables that are not columns of ",
            "data: ", paste0("'", unknown_targets, "'", collapse = ", "),
            ". Exclusions for them will have no effect.")
  }

  unknown_vars = setdiff(unique(unlist(adjustment_exclusions)), data_names)
  if (length(unknown_vars) > 0L) {
    warning("adjustment_exclusions asks to exclude variables that are not ",
            "columns of data: ",
            paste0("'", unknown_vars, "'", collapse = ", "),
            ". They will not be excluded from any adjustment set.")
  }

  invisible(NULL)
}
