# Helpers for the adjustment_exclusions argument of varimpact().
#
# varimpact() expands the analysis data before it becomes an adjustment matrix,
# so one input variable can appear in W under several different names:
#
#   a numeric V1                  -> "V1"
#   a factor  F1                  -> one dummy per level: "F1XXb", "F1XXc"
#                                    (see factors_to_indicators())
#   either one, if it had NAs     -> "Imiss_V1"
#
# Matching an exclusion request against colnames(W) literally would therefore
# remove a numeric variable but silently leave every dummy of a factor in place.

# Return the columns of an adjustment matrix that belong to one input variable.
adjustment_columns = function(name, w_names) {
  w_names[w_names == name |
          startsWith(w_names, paste0(name, "XX")) |
          w_names == paste0("Imiss_", name)]
}

# Drop the adjustment variables the user excluded for this candidate variable.
# Returns W unchanged when nothing was requested for it.
#
# This runs inside the per-variable, per-fold loop, so it stays silent: the
# spelling of adjustment_exclusions is checked once up front by
# check_adjustment_exclusions() instead of warning on every fold.
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

# Validate adjustment_exclusions once, against the columns of the analysis data.
#
# A name that matches nothing is almost always a typo, and excluding nothing
# silently is the dangerous failure here: the user would believe they had
# adjusted away a mediator or a collider when they had not. Warn rather than
# stop, so that a stale entry left over from an earlier version of the data
# does not abort a long run.
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
