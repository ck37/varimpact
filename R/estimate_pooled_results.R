# Qbounds: the bounds that were used to map the outcome into [0, 1] before the
# CV-TMLE was run, i.e. the same vector varimpact() passed down to
# apply_tmle_to_validation(). For a binary outcome this is c(0, 1) and every
# transformation below is the identity. For a continuous outcome it is the
# (10%-widened) range of Y, and we use it to map the results back onto the
# scale of the original outcome. See issue #8.
estimate_pooled_results = function(fold_results,
                                   fluctuation = "logistic",
                                   verbose = FALSE,
                                   Qbounds = c(0, 1)) {
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
    # A fold that produced no validation predictions for this bin leaves
    # val_preds as a zero-row dataframe rather than NULL - the fold ran, this
    # bin just had nothing in it. cbind()ing fold_num onto zero rows fails with
    # "arguments imply differing number of rows: 0, 1", so treat it as a failed
    # fold, which is what it is.
    if (is.null(fold$val_preds) || NROW(fold$val_preds) == 0L) {
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

  # delta = 0 marks validation observations that are missing Y or A. Their
  # clever covariate (HAW) is already 0, so they contribute nothing to the
  # fluctuation or to the influence curve.
  if (is.null(data$delta)) {
    data$delta = 1
  }

  if (min(data$Q_hat) < 0 || max(data$Q_hat) > 1) {
    stop("estimate_pooled_results(): predicted Q_hat values must lie in [0, 1]; ",
         "observed range ", paste(signif(range(data$Q_hat), 4), collapse = " to "), ".")
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
    if (inherits(data$logit_Q_hat, "try-error")) {
      stop("estimate_pooled_results(): qlogis() failed on Q_hat: ",
           conditionMessage(attr(data$logit_Q_hat, "condition")))
    }
    #}

    # Estimate epsilon
    if (verbose) cat("Estimating epsilon: ")

    if (fluctuation == "logistic") {
      suppressWarnings({
        #epsilon = coef(glm(Y_star ~ -1 + offset(logit_Q_hat) + H1W,
        #epsilon = coef(glm(Y_star ~ -1 + offset(logit_Q_hat) + HAW,
        #                 data = data, family = "binomial"))
        # offset() has to be called unqualified: stats::offset() is not
        # recognized as the formula's offset special, so logit_Q_hat would be
        # fit as an ordinary covariate and epsilon would come back with two
        # elements instead of one.
        reg = try(stats::glm(Y_star ~ -1 + offset(logit_Q_hat) + HAW,
                  data = data, family = "binomial",
                  subset = data$delta == 1))
        if ("try-error" %in% class(reg)) {
          stop("estimate_pooled_results(): the fluctuation regression for epsilon ",
               "failed: ", conditionMessage(attr(reg, "condition")))
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
      stop("estimate_pooled_results(): could not extract epsilon from the ",
           "fluctuation regression: ", conditionMessage(attr(epsilon, "condition")))
    } else {

      if (verbose) cat("Fluctuating Q_star\n")

      # Fluctuate Q to get Q_star
      Q_star = data$logit_Q_hat + epsilon * data$H1W
      #Q_star = data$logit_Q_hat + epsilon * data$HAW

      if (verbose) cat("Transforming Q_star\n")
      #if (length(unique(data$Y)) == 2) {
      Q_star = plogis(Q_star)
      #}

      # Map back onto the scale of the original outcome.
      #
      # apply_tmle_to_validation() ran the whole CV-TMLE on
      # Y_star = (Y - Qbounds[1]) / diff(Qbounds), so both Q_star and Y_star
      # currently live on that [0, 1] scale. The treatment-specific mean is
      # linear in Y, so its inverse is just the inverse of that map, and every
      # quantity built from Q_star and Y_star below inherits the right scale:
      #
      #   theta_original = theta_star * diff(Qbounds) + Qbounds[1]
      #   IC_original    = diff(Qbounds) * IC_star
      #
      # The location shift cancels out of the influence curve, because both of
      # its terms are differences: HAW * (Y_star - Q_star), and
      # Q_star - mean(Q_star).
      #
      # For a binary outcome Qbounds is c(0, 1) and this is a no-op, which is
      # why it is applied unconditionally rather than behind a family check.
      if (!is.null(Qbounds) && length(Qbounds) == 2L) {
        if (verbose && !identical(as.numeric(Qbounds), c(0, 1))) {
          cat("Mapping Q_star back to the outcome scale using Qbounds:",
              Qbounds, "\n")
        }
        Q_star = Q_star * diff(Qbounds) + Qbounds[1]
        # Observations with delta == 0 carry a placeholder Y_star of 0; they get
        # zero weight through HAW, so their rescaled value is irrelevant.
        data$Y_star = data$Y_star * diff(Qbounds) + Qbounds[1]
      }

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
        #with(fold_data, (A / g1W_hat) * (Y - Q_star) + Q_star - theta)
        # HAW = A * delta / (g1W * g.Delta), so observations missing Y or A
        # drop out of the residual term rather than turning it into an NA.
        result = with(fold_data, HAW * (Y_star - Q_star) +
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
          cat("gAW_total zeros:", sum(data$gAW_total == 0), "\n")
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

  # tapply() and by() key on the fold numbers that actually appear in the data,
  # so a fold contributing no rows for this bin gets no slot - and every later
  # fold shifts down a position. Callers index these BY FOLD NUMBER
  # (thetas[fold], influence_curves[[fold]]), so that shift either runs off the
  # end or, worse, silently attributes one fold's estimate to another. Re-expand
  # to one slot per fold, empty where the fold contributed nothing.
  num_folds = length(fold_results)
  align_by_fold = function(x, empty) {
    present = suppressWarnings(as.integer(names(x)))
    keep = !is.na(present) & present >= 1L & present <= num_folds
    out = rep(empty, num_folds)
    out[present[keep]] = x[keep]
    out
  }
  if (!is.null(thetas)) {
    # as.numeric() drops the fold-number names that align_by_fold() reads, so
    # put them back.
    thetas_flat = as.numeric(thetas)
    names(thetas_flat) = names(thetas)
    thetas = align_by_fold(thetas_flat, NA_real_)
  }
  if (!is.null(influence_curves)) {
    # as.list() strips the "by" class while keeping the fold-number names.
    influence_curves = align_by_fold(as.list(influence_curves), list(NULL))
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
