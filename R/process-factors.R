process_factors = function(data.fac,
                           quantile_probs_factor,
                           miss.cut,
                           verbose = FALSE) {

  #####################
  if (ncol(data.fac) > 0L) {


    if (verbose) cat("Processing factors. Start count:", ncol(data.fac), "\n")

    ######################################
    # Replace blank factor values with NA's.

    # We re-use this num_cols variable in the next section.
    num_cols = ncol(data.fac)
    for (i in 1:num_cols) {
      new_factor = as.character(data.fac[, i])
      # The exclude argument replaces any empty strings with NAs.
      new_factor = factor(new_factor, exclude = "")
      data.fac[, i] = new_factor
    }

    ###################
    # For each factor, apply function and get rid of those where
    # 'true' data.fac is data frame of variables that are factors
    if (!is.null(quantile_probs_factor)) {
      data.fac = restrict_by_quantiles(data.fac, quantile_probs = quantile_probs_factor)
    }

    dropped_cols = num_cols - ncol(data.fac)

    if (verbose) {
      if (dropped_cols > 0) {
        cat("Dropped", dropped_cols, "factors due to lack of variation.\n")
      } else {
        cat("No factors dropped due to lack of variation.\n")
      }
    }

    # We don't seem to use this yet.
    # num.cat = sapply(data.fac, length_unique)

    ######################
    # Remove columns with missing data % greater than the threshold.
    sum_nas = sapply(data.fac, sum_na)

    if (verbose) cat("Factors with missingness:", sum(sum_nas > 0L), "\n")

    miss_pct = sum_nas / nrow(data.fac)

    data.fac = data.fac[, miss_pct < miss.cut, drop = FALSE]

    if (verbose) {
      cat("Dropped", sum(miss_pct >= miss.cut), "factors due to the missingness threshold.\n")
    }

    # Save how many separate factors we have in this dataframe.
    num_factors = ncol(data.fac)

    factor_results = factors_to_indicators(data.fac, verbose = verbose)

    datafac.dum = factor_results$data
    # Here 1 = defined, 0 = missing.
    miss.fac = factor_results$missing_indicators

    if (verbose) {
      cat("End factor count:", num_factors, "Indicators:", ncol(datafac.dum),
          "Missing indicators:", ncol(miss.fac), "\n")
    }
  } else {
    num_factors = 0L
    miss.fac = NULL
    datafac.dum = NULL
    data.fac = NULL
  }

  (results =
    list(
      num_factors = num_factors,
      miss.fac = miss.fac,
      datafac.dum = datafac.dum,
      data.fac = data.fac
  ))

}
