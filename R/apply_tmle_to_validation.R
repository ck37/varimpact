apply_tmle_to_validation = function(Y,
                                    A,
                                    W,
                                    family,
                                    delta = rep(1, length(Y)),
                                    tmle,
                                    id = 1:length(Y),
                                    verbose = F) {

  ###########
  # CK: actually, I think we don't even need to use Y
  # Transform Y to Y_star
  #Y_star = Y
  #if (tmle$map_to_ystar) {
  #  Y_star = (Y_star - tmle$ab[1]) / diff(tmle$ab)
  #}

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

  # TODO: return Y_star rather than Y, for use in fluctuation step?
  data = data.frame(Y = Y, A = A, Q_hat = Q_hat,
                    g1W_hat = g1W_hat,
                    H1W = H1W)

  results = data

  return(results)
}
