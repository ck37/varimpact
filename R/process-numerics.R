process_numerics =
  function(data.num,
           quantile_probs_numeric,
           miss.cut,
           bins_numeric,
           impute,
           verbose = FALSE) {

  n = nrow(data.num)

  if (ncol(data.num) > 0) {
    num_cols = ncol(data.num)
    if (verbose) cat("Processing numerics. Start count:", num_cols, "\n")

    # Remove columns where the 0.1 and 0.9 quantiles have the same value, i.e. insufficent variation.
    # TODO: set this is a configurable setting?
    if (!is.null(quantile_probs_numeric)) {
      data.num = restrict_by_quantiles(data.num, quantile_probs = quantile_probs_numeric)
    }

    if (verbose) {
      num_dropped = num_cols - ncol(data.num)
      if (num_dropped > 0) {
        cat("Dropped", num_dropped, "numerics due to lack of variation.\n")
      } else {
        cat("No numerics dropped due to lack of variation.\n")
      }
    }

    # Save how many numeric variables we have in this dataframe.
    num_numeric = ncol(data.num)
  } else {
    num_numeric = 0L
  }

  if (num_numeric > 0) {
    if (verbose) cat("Cleaning up", num_numeric, "numeric variables.\n")
    # Make deciles for continuous variables
    X = data.num
    xc = dim(X)[2]
    # We don't seem to use this qt variable.
    qt = apply(na.omit(X), 2, quantile, probs = seq(0.1, 0.9, 0.1))
    newX = NULL
    coln = NULL
    varn = colnames(X)

    num.cat = apply(X, 2, length_unique)

    # Matrix to store the numeric columns converted to bin levels (integer value per quantile).
    numerics_binned = matrix(nrow = n, ncol = num_numeric)

    # Make a list to store the levels for each numeric variable.
    numeric_levels = vector("list", num_numeric)

    for (numeric_i in 1:num_numeric) {
      # Because we do not specify "drop" within the brackets, Xt is now a vector.
      Xt = X[, numeric_i]

      name = colnames(data.num)[numeric_i]

      if (verbose) {
        cat("Processing", name, numeric_i, "of", num_numeric, "\n")
      }

      # Suppress the warning that can occur when there are fewer than the desired
      # maximum number of bins, as specified by bins_numeric. We should be able to
      # see this as var_binned containing fewer than bins_numeric columns.
      # Warning is in .cut2(): min(xx[xx > upper])
      # "no non-missing arguments to min; returning Inf"
      num_unique_vals = length(setdiff(unique(Xt), NA))

      num_breaks = min(bins_numeric, num_unique_vals)

      arules_method = "frequency"

      # No need to apply tiling, we already have a limited # of unique vals.
      # TODO: return info on this, and possible a notice message.
      if (num_unique_vals <= bins_numeric) {
        arules_method = "interval"
      }

      suppressWarnings({
        # Discretize into up to 10 quantiles (by default), configurable based on
        # bins_numeric argument.
        # This returns a factor version of the discretized variable.
        tryCatch({ var_binned_names = arules::discretize(Xt,
                                              method = arules_method,
                                              breaks = num_breaks,
                                              ordered = TRUE)
        }, error = function(error) {
          # This can happen with skewed distributions where multiple breaks are not unique.
          print(error)
          cat("Error: could not discretize numeric", numeric_i, "", name, "\n")
          cat("Unique values:", length(unique(Xt)), "\n")
          cat("Switching to cluster-based discretization.\n")
          tryCatch({
          var_binned_names = arules::discretize(Xt, method = "cluster",
                                                breaks = num_breaks,
                                                ordered = TRUE)},
            error = function(error2) {
              # TODO: use another package/function to discretize.
              print(error2)
              cat("Cluster-based discretization failed - using all levels.")
              var_binned_names = factor(Xt)
            })

        })
      })

      # Save the levels for future usage.
      numeric_levels[[numeric_i]] = levels(var_binned_names)
      # This converts the factor variable to just the quantile numbers.
      var_binned = as.numeric(var_binned_names)
      numerics_binned[, numeric_i] =  var_binned

      if (verbose) {
        cat("New levels:", paste(levels(var_binned_names), collapse = ", "), "\n")
      }
    }
    colnames(numerics_binned) = varn
    data.cont.dist = numerics_binned

    ###############
    # Missing Basis for numeric variables, post-binning.

    n.cont = nrow(data.cont.dist)

    sum_nas = apply(data.cont.dist, 2, sum_na)
    nmesX = colnames(data.cont.dist)
    miss.cont = NULL
    nmesm = NULL

    # Create imputed version of the numeric dataframe.
    # This is used as the adjustment set, but not used when generating the treatment assignment vector.
    data.numW = data.num

    # Loop over each binned numeric variable.
    # TODO: do this as part of the binning process.
    for (k in 1:num_numeric) {
      # Check if that variable has any missing values.
      if (sum_nas[k] > 0) {
        # The effect is that the basis is set to 1 if it exists and 0 if it's missing.
        ix = as.numeric(!is.na(data.cont.dist[, k]))
        miss.cont = cbind(miss.cont, ix)
        # TODO: convert to paste0
        nmesm = c(nmesm, paste("Imiss_", nmesX[k], sep = ""))
      }
    }
    # if(is.null(miss.cont)){miss.cont= rep(1,n.cont)}
    colnames(miss.cont) = nmesm

    # Impute missing data in numeric columns.
    if (impute == "zero") {
      data.numW[is.na(data.num)] = 0
      impute_info = 0
    } else if (impute == "median") {
      impute_info = caret::preProcess(data.num, method = "medianImpute")
      data.numW = predict(impute_info, data.num)
    } else if (impute == "mean") {
      stop("Mean imputation not implemented yet. Please use another imputation method.")
    } else if (impute == "knn") {
      # NOTE: this also results in caret centering and scaling the data.
      impute_info = caret::preProcess(data.num, method = "knnImpute")
      data.numW = predict(impute_info, data.num)
    }

    # Confirm that there are no missing values remaining in data.numW
    stopifnot(sum(is.na(data.numW)) == 0)
  } else {
    data.cont.dist = NULL
    miss.cont = NULL
    data.numW = NULL
    impute_info = NULL
    data.num = NULL
  }

  (results =
      list(
        data.cont.dist = data.cont.dist,
        num_numeric = num_numeric,
        miss.cont = miss.cont,
        data.num = data.num,
        data.numW = data.numW,
        impute_info = impute_info
  ))
}
