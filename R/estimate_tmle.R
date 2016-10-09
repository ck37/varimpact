#' Get TMLE estimate.
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
  if (is.null(delta) == F) {
    ss = sum(delta == 0)
    if (ss > 0 & ss < 10) {
      inc[delta == 0] = FALSE
    }
  }
  Y = Y[inc]
  A = A[inc]
  W = W[inc, , drop = F]
  delta = delta[inc]
  tmle.1 = tmle::tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                      Q.SL.library = Q.lib, family = family, verbose = verbose)
  g1 = tmle.1$g$g1W
  Qst = tmle.1$Qstar[, 2]
  theta = mean(Qst)
  # Influence curve
  IC = (A / g1) * (Y - Qst) + Qst - theta
  return(list(theta = theta, IC = IC))
}
