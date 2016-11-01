estimate_pooled_results = function(fold_results, fluctuation = "logistic", verbose = F) {
  # Fold results is a list with results from each fold.

  # Each fold result should have at least this element:
  # val_preds dataframe, with Y_star, g, Q, H.

  data = do.call(rbind, lapply(1:length(fold_results), function(i) {
    fold = fold_results[[i]]
    # Save the fold number so we can use it to generate fold-specific estimates.
    if (is.null(fold$val_preds)) {
      NULL
    } else {
      df = cbind(fold$val_preds, fold_num = i)
      df
    }
  }))


  if (!is.null(data)) {
    n = nrow(data)

    # If Y is binary, take logit of Q.
    #if (length(unique(data$Y)) == 2) {
    data$Q_hat = qlogis(data$Q_hat)
    #}

    # Estimate epsilon
    if (fluctuation == "logistic") {
      #epsilon = coef(glm(Y_star ~ -1 + offset(Q_hat) + H1W, data = data, family = "binomial"))
      # Use more stable version where clever covariate is the weight, and now we
      # have an intercept. Causal 2, Lecture 3, slide 51.
      # We have to suppressWarnings about "non-integrate #successes in binomial glm".
      suppressWarnings({
        epsilon = coef(glm(Y_star ~ offset(Q_hat), weights = H1W, data = data, family = "binomial"))
      })
    } else {
      stop("Only support logistic fluctuation currently.")
      # TBD.
    }

    # Fluctuate Q to get Q_star
    Q_star = data$Q_hat + epsilon * data$H1W

    #if (length(unique(data$Y)) == 2) {
    Q_star = plogis(Q_star)
    #}

    # Estimate parameter on every validation fold.
    thetas = tapply(Q_star, data$fold_num, mean, na.rm = T)

    # Take average across folds to get final estimate.
    #theta = mean(thetas)

    # Get influence curve per fold.
    # Influence_curves here is a matrix, with one column per fold.
    influence_curves = simplify2array(base::by(data, data$fold_num, function(fold_data) {
      #with(fold_data, (A / g1W_hat) * (Y - Q_star) + Q_star - theta)
      with(fold_data, (A / g1W_hat) * (Y_star - Q_star) + Q_star - mean(Q_star, na.rm=T))
    }))

    # Old version:
    #influence_curve = with(data, (A / g1W_hat) * (Y - Q_star) + Q_star - theta)

    # Calculate standard error.
    #std_err = stats::var(influence_curves) / n
  } else {
    # All folds must have failed.
    thetas = NULL
    influence_curves = NULL
    epsilon = NULL

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
