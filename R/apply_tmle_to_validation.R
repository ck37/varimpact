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
    Y_star = (Y_star - tmle$Qbounds[1]) / diff(tmle$Qbounds)
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

  # Predict Q
  sl_pred = predict(tmle$q_model, data.frame(Z, A, W), onlySL = T)
  Q_hat = sl_pred$pred

  # Predict g
  sl_pred = predict(tmle$g_model, W, type = "response", onlySL = T)
  g1W_hat = sl_pred$pred

  # Bound g1W_hat
  # TODO: double-check this with Alan.
  g1W_hat = .bound(g1W_hat, tmle$gbounds)

  # Create clever covariate.
  H1W = (A == 1) / g1W_hat

  # TODO: handle delta in some way?

  # Return results

  # We return Y_star rather than Y, for use in pooled fluctuation step.
  data = data.frame(Y_star = Y_star, A = A, Q_hat = Q_hat,
                    g1W_hat = g1W_hat,
                    H1W = H1W)

  results = data

  return(results)
}
