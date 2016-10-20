#' Get TMLE estimate: E[Y | A = 1, W].
#'
#' @param Y Outcome variable
#' @param A Treatment indicator
#' @param W Dataframe of adjustment covariates
#' @param family Binomial or gaussian
#' @param delta Indicator of missing outcome or treatment assignment. 1 - observed, 0 - missing.
#' @param Q.lib SuperLearner library for estimating Q (potential outcome)
#' @param g.lib SuperLearner library for estimating g (propensity score)
#' @param verbose If true output extra information during execution.
#' @importFrom tmle tmle
estimate_tmle = function(Y,
                         A,
                         W,
                         family,
                         delta = NULL,
                         Q.lib,
                         g.lib,
                         verbose = F) {

  if (!family %in% c("binomial", "gaussian")) {
    stop('Estimate_tmle: family must be either "binomial" or "gaussian".')
  }

  # Because of quirk of program, delete observations with delta=0 if #>0
  # & < 10
  n = length(Y)
  inc = rep(TRUE, n)

  if (!is.null(delta)) {
    ss = sum(delta == 0)
    if (ss > 0 & ss < 10) {
      inc[delta == 0] = FALSE
    }
  }

  Y = Y[inc]
  A = A[inc]
  W = W[inc, , drop = F]

  delta = delta[inc]

  # Here we are using tmle but not using the treatment effect estimate.
  # We're actually using the underlying variables to estimate Y_a.
  tmle.1 = tmle::tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                      Q.SL.library = Q.lib, family = family, verbose = verbose)

  # Propensity score for treatment.
  g1 = tmle.1$g$g1W

  # Unit's estimated outcome under treatment: hat(Y) | A = 1, W
  # This is after the fluctuation step, so it is targeted.
  Qst = tmle.1$Qstar[, 2]

  # E[Y | A = 1, W]
  theta = mean(Qst)

  # Influence curve
  IC = (A / g1) * (Y - Qst) + Qst - theta

  # Compile results.
  result = list(theta = theta, IC = IC)

  return(result)
}
