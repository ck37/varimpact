apply_tmle_to_validation = function(Y,
                                    A,
                                    W,
                                    family,
                                    delta = rep(1, length(Y)),
                                    tmle,
                                    id = 1:length(Y),
                                    verbose = F) {

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
      cat("New Y_star range:", range(Y_star), "\n")
    }
  }

  if (max(Y_star) > 1 | min(Y_star) < 0) {
    cat("Error on Y_star's range for logistic fluctuation\n")
    cat("Y_star distribution:\n")
    print(summary(Y_star))
    cat("Qbounds:", tmle$Qbounds, "\n")
    cat("Stage1 Qbounds:", tmle$stage1_Qbounds, "\n")
    cat("Stage1 ab:", tmle$ab, "\n")
    cat("Validation Y range:", range(Y), "\n")
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
    browser()
    stop("apply_tmle_to_validation() failed during prediction of Q(1, W).")
  })

  if (verbose) cat("Bounding Q_hat to", tmle$stage1_Qbounds, "\n")
  Q_hat = .bound(Q_hat, tmle$stage1_Qbounds)

  if (min(Q_hat) < 0 || max(Q_hat) > 1) {
    cat("Error: predicted Q_hat outside of [0, 1] bounds.\n")
    browser()
  }

  # Predict g
  tryCatch({
    sl_pred = predict(tmle$g_model, W, type = "response", onlySL = T)
    g1W_hat = sl_pred$pred
  }, error = function(e) {
    print(e)
    print(tmle$g_model)
    browser()
    stop("apply_tmle_to_validation() failed during prediction of g.")
  })

  if (verbose) cat("Current range of g1W on test:", range(g1W_hat), "\n")

  # Truncate g1W_hat
  # TODO: double-check this with Alan.
  g1W_hat_truncated = .bound(g1W_hat, tmle$gbounds)
  if (verbose) cat("Truncating g1W on test using bounds:", tmle$gbounds, "\n")

  # Create clever covariate.
  # H1W = (A == 1) / g1W_hat_truncated
  # Based on Jeremy Coyle's CV-TMLE implementation.
  H1W = 1 / g1W_hat_truncated
  HAW = A * H1W

  # TODO: handle delta in some way?

  if (verbose) {
    cat("Mean Y_a on validation:", mean(Y[A == 1]), "\n")
    cat("Mean Y_star_a on validation:", mean(Y_star[A == 1]), "\n")
    cat("Mean Q_bar_a on validation:", mean(Q_hat[A == 1]), "\n")
  }

  ####################
  # Return results

  # We return Y_star rather than Y, for use in pooled fluctuation step.
  data = data.frame(Y_star = Y_star, A = A, Q_hat = Q_hat,
                    g1W_hat = g1W_hat_truncated,
                    H1W = H1W,
                    HAW = HAW)

  results = data

  return(results)
}
