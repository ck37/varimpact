estimate_pooled_results = function(fold_results,
                                   fluctuation = "logistic",
                                   verbose = FALSE) {
  # Fold results is a list with test results from each fold.

  # Each fold result should have at least this element:
  # val_preds dataframe, with Y_star, g, Q, H.
  # TODO: need to change this check, it doesn't work correctly.
  #num_fails = sum(is.null(sapply(fold_results, `[[`, "level")))
  #if (verbose) {
  #  cat("Number of fold failures:", num_fails, "of", length(fold_results), "\n")
  #}

  # Placeholder results to return in case of error.
  results = list(
    thetas = NULL,
    influence_curves = NULL,
    epsilon = NULL
  )

  #if (num_fails == length(fold_results)) {
    # Every fold failed.
  #  if (verbose) cat("Error: every fold failed.\n")
  #  return(results)
  #}

  # browser()

  # Extract the results from each CV-TMLE fold and rbind into a single dataframe.
  data = do.call(rbind, lapply(1:length(fold_results), function(i) {
    fold = fold_results[[i]]
    # Save the fold number so we can use it to generate fold-specific estimates.
    if (is.null(fold$val_preds)) {
      # Skip folds that failed.
      NULL
    } else {
      # val_preds is a dataframe with columns: Y_star, g, Q, H
      df = cbind(fold$val_preds, fold_num = i)
      df
    }
  }))

  if (is.null(data)) {
    # Every fold failed.
    if (verbose) cat("Error: every fold failed.\n")
    return(results)
  }

  if (min(data$Q_hat) < 0 || max(data$Q_hat) > 1) {
    cat("Error: some predicted values of Q_hat are out of bounds.",
        "They should be in [0, 1].\n")
    print(summary(data$Q_hat))
    browser()
  }

  # Set some default values in case of a future error.
  thetas = NULL
  influence_curves = NULL
  epsilon = NULL


  if (!is.null(data)) {
    n = nrow(data)

    # If Y is binary, take logit of Q.
    #if (length(unique(data$Y)) == 2) {

    # Look at thetas prior to fluctuation.
    pre_thetas = tapply(data$Q_hat, data$fold_num, mean, na.rm = TRUE)
    if (verbose) cat("Pre-fluctuation thetas:", pre_thetas, "\n")

    # If Q is binary or continuous we still want to take logit of predicted values.
    # See tmle::estimateQ where it does this after predicting Q.
    data$logit_Q_hat = try(stats::qlogis(data$Q_hat))
    if (class(data$logit_Q_hat) == "try-error") {
      cat("Error in estimate_pooled_results() with qlogis()\n")
      print(summary(data$Q_hat))
      browser()
    }
    #}

    # Estimate epsilon
    if (verbose) cat("Estimating epsilon: ")

    if (fluctuation == "logistic") {
      suppressWarnings({
        #epsilon = coef(glm(Y_star ~ -1 + offset(logit_Q_hat) + H1W,
        #epsilon = coef(glm(Y_star ~ -1 + offset(logit_Q_hat) + HAW,
        #                 data = data, family = "binomial"))
        reg = try(stats::glm(Y_star ~ -1 + stats::offset(logit_Q_hat) + HAW,
                  data = data, family = "binomial"))
        if ("try-error" %in% class(reg)) {
          cat("Error in epsilon regression.\n")
          browser()
        }
        epsilon = try(stats::coef(reg))
      })
      # Use more stable version where clever covariate is the weight, and now we
      # have an intercept. Causal 2, Lecture 3, slide 51.
      # We have to suppressWarnings about "non-integrate #successes in binomial glm".
      #suppressWarnings({
        # Catch an error if one occurs here.
        #epsilon = try(coef(glm(Y_star ~ offset(logit_Q_hat),
      #  epsilon = try(coef(glm(Y_star ~ .,
      #                         offset = logit_Q_hat,
      #                         weights = H1W,
      #                         data = data, family = "binomial")))
      #})
      if (verbose) cat(epsilon, "\n")
    } else {
      # No need to support linear fluctuation as it does not respect model bounds.
      stop("Only support logistic fluctuation currently.")
      # TBD.
    }

    if ("try-error" %in% class(epsilon)) {
      if (verbose) cat("Error when estimating epsilon.\n")
      print(summary(data$Y_star))
      browser()
    } else {

      if (verbose) cat("Fluctuating Q_star\n")

      # Fluctuate Q to get Q_star
      Q_star = data$logit_Q_hat + epsilon * data$H1W
      #Q_star = data$logit_Q_hat + epsilon * data$HAW

      if (verbose) cat("Transforming Q_star\n")
      #if (length(unique(data$Y)) == 2) {
      Q_star = plogis(Q_star)
      #}

      if (verbose) cat("Estimating per-fold thetas: ")

      # Estimate treatment-specific mean parameter on every validation fold.
      thetas = tapply(Q_star, data$fold_num, mean, na.rm = TRUE)
      if (verbose) cat(thetas, "\n")

      # Take average across folds to get final estimate.
      #theta = mean(thetas)

      # Move Q_star into the data so that it can be analyzed per-fold.
      data$Q_star = Q_star
      rm(Q_star)

      if (verbose) cat("Calculating per-fold influence curves\n")

      # Get influence curve per fold - for treatment-specific mean.
      # Influence_curves here is a list, where each element is a result.
      # We can't convert to a matrix because lengths are different.
      # TODO: figure out why this can generate NaNs
      influence_curves = base::by(data, data$fold_num, function(fold_data) {
        if (F && verbose) {
          with(fold_data,
               cat("A:", length(A), "g1W_hat:", length(g1W_hat), "Y_star:", length(Y_star),
                   "Q_star:", length(Q_star), "\n"))
        }
        #with(fold_data, (A / g1W_hat) * (Y - Q_star) + Q_star - theta)
        result = with(fold_data, (A / g1W_hat) * (Y_star - Q_star) +
                        Q_star - mean(Q_star, na.rm = TRUE))
        #if (verbose) cat("Result:", class(result), "Length:", length(result), "\n")
        result
      })

      # Check for NaNs.
      num_nans = sum(sapply(influence_curves, function(curve) sum(is.nan(curve))))
      if (num_nans > 0) {
        if (verbose) {
          cat("Error: influence curves contain", num_nans, "NaNs.\n")
          cat("g1W_hat zeros:", sum(data$g1W_hat == 0), "\n")
        }
      }

      #if (verbose) cat("IC class:", class(influence_curves), "\n")

      # Old version:
      #influence_curve = with(data, (A / g1W_hat) * (Y - Q_star) + Q_star - theta)

      # Calculate standard error.
      #std_err = stats::var(influence_curves) / n
    }
  }

  if (is.null(thetas))  {
    # All folds must have failed.
    if (verbose) cat("No pooled results. All folds seemed to have failed.\n")
  }

  # Compile results
  results = list(
    #theta = theta,
    thetas = thetas,
    influence_curves = influence_curves,
    #std_err = std_err,
    epsilon = epsilon
  )

  return(results)
}
