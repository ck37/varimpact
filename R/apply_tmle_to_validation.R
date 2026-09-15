# Predict the missingness mechanism P(Delta = 1 | Z = 0, A = 1, W) on new data,
# using the model that estimate_tmle2() fit on the training fold.
# Returns a vector of 1s if the training fold had no missingness, in which case
# no model was fit.
predict_g_delta = function(model, W) {
  if (is.null(model)) {
    # No missingness in the training fold, so P(Delta = 1 | A, W) = 1.
    return(rep(1, nrow(W)))
  }

  # The missingness model was fit on data.frame(Delta, Z, A, W), so the
  # prediction data needs those same columns. We want the counterfactual in
  # which the observation is treated, matching the "Z0A1" column that
  # estimate_tmle2() uses on the training data.
  newdata = data.frame(Z = 0, A = 1, W)

  # Z is constant in the data the model was fit on, so the fit can be
  # rank-deficient and predict() warns. estimate_tmle2() suppresses the same
  # warning when it predicts this model on the training data.
  tryCatch(suppressWarnings({
    if (inherits(model, "SuperLearner")) {
      as.vector(predict(model, newdata, onlySL = TRUE)$pred)
    } else {
      as.vector(stats::predict(model, newdata = newdata, type = "response"))
    }
  }), error = function(e) {
    print(e)
    stop("apply_tmle_to_validation() failed during prediction of g.Delta.")
  })
}

apply_tmle_to_validation =
  function(Y,
           A,
           W,
           family,
           delta = rep(1, length(Y)),
           tmle,
           id = 1:length(Y),
           verbose = FALSE) {

  # delta = 1 when both Y and A are observed, 0 otherwise. Observations with
  # delta = 0 contribute no residual to the fluctuation step or the influence
  # curve; they are instead accounted for by the missingness mechanism
  # g.Delta that was estimated on the training fold.
  delta = as.numeric(delta)
  observed = delta == 1

  if (sum(observed) == 0L) {
    stop("apply_tmle_to_validation(): no observations have both Y and A observed.")
  }

  ###########
  # Transform Y to Y_star, needed for later fluctuation step.
  Y_star = Y
  if (tmle$map_to_ystar) {
    # Use the Qbounds from the full range of Y,
    # not the tmle$ab that is based only on training data.
    # Y_star = (Y_star - tmle$ab[1]) / diff(tmle$ab)
    # If Y is binary this will have no effect as bounds & outcome are already {0, 1}.
    Y_star = (Y_star - tmle$Qbounds[1]) / diff(tmle$Qbounds)

    if (verbose) {
      cat("Mapped Y to Y_star using Qbounds:", tmle$Qbounds, "\n")
      cat("New Y_star range:", range(Y_star[observed]), "\n")
    }
  }

  # Only the observed outcomes need to respect the bounds; the missing ones are
  # given zero weight below.
  if (any(is.na(Y_star[observed]))) {
    cat("Error: Y is missing for observations that are flagged as observed.\n")
    stop("delta must be 0 for any observation with a missing Y or A.")
  }

  if (max(Y_star[observed]) > 1 | min(Y_star[observed]) < 0) {
    cat("Error on Y_star's range for logistic fluctuation\n")
    cat("Y_star distribution:\n")
    print(summary(Y_star[observed]))
    cat("Qbounds:", tmle$Qbounds, "\n")
    cat("Stage1 Qbounds:", tmle$stage1_Qbounds, "\n")
    cat("Stage1 ab:", tmle$ab, "\n")
    cat("Validation Y range:", range(Y[observed]), "\n")
    stop("Y values must be 0 <= y <= 1")
  }

  # We only include this because TMLE functions use Z.
  Z = rep(0, length(Y))

  q_df = data.frame(Z, A = 1, W)

  # Predict Q(1, W)
  tryCatch({
    sl_pred = predict(tmle$q_model, q_df, onlySL = T)
    Q_hat = sl_pred$pred
  }, error = function(e) {
    print(e)
    print(tmle$q_model)
    #browser()
    stop("apply_tmle_to_validation() failed during prediction of Q(1, W).")
  })

  if (verbose) cat("Bounding Q_hat to", tmle$stage1_Qbounds, "\n")
  Q_hat = .bound(Q_hat, tmle$stage1_Qbounds)

  if (min(Q_hat) < 0 || max(Q_hat) > 1) {
    cat("Error: predicted Q_hat outside of [0, 1] bounds.\n")
    #browser()
  }

  # Predict g
  tryCatch({
    # Check specifically for a g_model that doesn't exist.
    if (is.null(tmle$g_model)) {
      stop("tmle$g_model has class = NULL")
    }
    sl_pred = predict(tmle$g_model, W, type = "response", onlySL = TRUE)
    g1W_hat = sl_pred$pred
  }, error = function(e) {
    print(e)
    print(tmle$g_model)
    #browser()
    stop("apply_tmle_to_validation() failed during prediction of g.")
  })

  if (verbose) cat("Current range of g1W on test:", range(g1W_hat), "\n")

  # Truncate g1W_hat
  # TODO: double-check this with Alan.
  g1W_hat_truncated = .bound(g1W_hat, tmle$gbounds)
  if (verbose) cat("Truncating g1W on test using bounds:", tmle$gbounds, "\n")

  # Predict the missingness mechanism, P(Delta = 1 | A = 1, W).
  gDelta_hat = predict_g_delta(tmle$g_delta_model, W)

  if (verbose) {
    cat("Range of g.Delta on test:", range(gDelta_hat), "\n")
    cat("Observations missing Y or A on test:", sum(!observed), "of", length(Y), "\n")
  }

  # Create clever covariate.
  # H1W = (A == 1) / g1W_hat_truncated
  # Based on Jeremy Coyle's CV-TMLE implementation.
  # The denominator includes the missingness mechanism so that observations
  # with an observed outcome are up-weighted to represent those without one.
  # We truncate the product rather than each term, matching estimate_tmle2().
  gAW_total = .bound(g1W_hat * gDelta_hat, tmle$gbounds)
  H1W = 1 / gAW_total
  HAW = A * delta * H1W

  # Observations missing Y (or A) receive zero weight via HAW, so replace their
  # NA outcome with a placeholder to keep it from propagating into the pooled
  # fluctuation and influence curve calculations.
  Y_star[!observed] = 0

  if (verbose) {
    cat("Mean Y_a on validation:", mean(Y[A == 1 & observed]), "\n")
    cat("Mean Y_star_a on validation:", mean(Y_star[A == 1 & observed]), "\n")
    cat("Mean Q_bar_a on validation:", mean(Q_hat[A == 1]), "\n")
  }

  ####################
  # Return results

  # We return Y_star rather than Y, for use in pooled fluctuation step.
  data = data.frame(Y_star = Y_star, A = A, Q_hat = Q_hat,
                    g1W_hat = g1W_hat_truncated,
                    gDelta_hat = gDelta_hat,
                    gAW_total = gAW_total,
                    delta = delta,
                    H1W = H1W,
                    HAW = HAW)

  results = data

  return(results)
}
