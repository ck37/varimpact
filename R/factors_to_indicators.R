#' Convert a dataframe of factor into separate indicators.
#'
#' More details
#'
#' @param factor_df Dataframe consisting only of factor variables.
#' @param miss_name_prefix Starting name for each missing indicator.
#' @param verbose Set to T for detailed output.
#'
#' @return A list of results.
#'
#' @export
factors_to_indicators = function(factor_df, miss_name_prefix = "Imiss_",
                                 verbose = F) {

  # TODO: confirm that factor_df contains only factor variables.

  sum_nas = apply(factor_df, 2, sum_na)

  if (verbose) cat("Factors with missingness:", sum(sum_nas > 0), "\n")

  facnames = names(factor_df)

  nam.fac = function(x, name) {
    num_chars = nchar(x)
    out = paste(name, substr(x, 2, num_chars), sep = "XX")
    # Remove spaces in variable names.
    out = gsub(" ", "", out)
    return(out)
  }

  ############
  # Missing Basis for Factors

  factor_names = colnames(factor_df)
  miss.fac = NULL
  names_miss = NULL

  # Calculate the number of missing values for each column.
  sum_nas = apply(factor_df, 2, sum_na)

  # Loop over each column in our new indicator dataframe.
  # TODO: use which() to only loop over indices with missing data.
  for (col_k in 1:ncol(factor_df)) {
    if (sum_nas[col_k] > 0) {
      # if (verbose) cat("Missing data in", factor_names[col_k], "\n")
      # Again, we are flagging non-missing as 1 and missing as 0 here.
      indicator_column = as.numeric(is.na(factor_df[, col_k]) == F)

      miss.fac = cbind(miss.fac, indicator_column)

      names_miss = c(names_miss, paste0(miss_name_prefix, factor_names[col_k]))
    }
  }
  colnames(miss.fac) = names_miss

  newX = NULL
  factor_names = NULL

  # TODO: remove this line?
  options(na.action = "na.pass")

  # Loop over each factor variable.
  for (i in 1:ncol(factor_df)) {
    # Note: we want to keep this variable name as x so that the names are short
    # from model.matrix. We rely on this variable name being only 1 character.
    x = factor_df[, i]
    # CK: looks like we are omitting the first level?
    if (T || sum_nas[i] == 0) {
      omit_levels = -1
    } else {
      # if there is missing data, also omit the last level (NA)
      omit_levels = -1 * c(1, length(levels(x)))
    }

    # Convert to a series of indicators.
    indicators = model.matrix(~ x - 1)[, omit_levels]
    names = colnames(indicators)

    # Any remaining missing data is set to 0.
    remaining_nas = sum(sapply(indicators, function(col) sum(is.na(col))))
    if (remaining_nas > 0) {
      if (verbose) cat("Replacing", remaining_nas, "remaining nas with 0s.\n")
      indicators[is.na(indicators)] = 0
    }

    # CK: why do we need this option?
    if (is.null(names)) {
      names2 = facnames[i]
    } else {
      # Clean up the names for each indicator.
      names2 = nam.fac(names, facnames[i])
    }

    # Accumulate the names.
    factor_names = c(factor_names, names2)
    # Append the next set of indicators to the resulting matrix.
    newX = cbind(newX, indicators)
  }

  colnames(newX) = factor_names

  ##################
  # Indexing vector for dummy basis back to original factors.
  cc = regexpr("XX", factor_names)
  ncc = nchar(factor_names)
  cc[cc < 0] = ncc[cc < 0] + 1

  # NOTE: we don't seem to actually use this in varImpact() at the moment.
  factor_index = substr(factor_names, 1, cc - 1)

  # Create a list to hold multi-variable results.
  results = list(data = newX, factor_index = factor_index, missing_indicators = miss.fac)

  results
}
