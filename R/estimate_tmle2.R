#' Get TMLE estimate: E[Y | A = 1, W].
#'
#' @param Y Outcome variable
#' @param A Treatment indicator
#' @param W Dataframe of adjustment covariates
#' @param family Outcome family - binomial or gaussian
#' @param delta Indicator of missing outcome or treatment assignment. 1 - observed, 0 - missing.
#' @param Q.lib SuperLearner library for estimating Q (potential outcome)
#' @param g.lib SuperLearner library for estimating g (propensity score)
#' @param id Optional subject-level identifier.
#' @param Qbounds Bounds on Q
#' @param gbound Bounds on G
#' @param alpha TBD, from TMLE package
#' @param fluctuation Only logistic is currently supported.
#' @param V Number of folds for SuperLearner
#' @param verbose If true output extra information during execution.
#' @importFrom tmle tmle
#' @importFrom stats as.formula binomial coef glm plogis poisson predict qlogis
#' @importFrom utils packageDescription
#' @export
estimate_tmle2 =
  function(Y,
           A,
           W,
           family,
           delta = rep(1, length(Y)),
           Q.lib,
           g.lib,
           id = 1:length(Y),
           Qbounds = NULL,
           gbound = 0.025,
           alpha = 0.995,
           fluctuation = "logistic",
           V = 10,
           verbose = F) {

  if (!family %in% c("binomial", "gaussian")) {
    stop('Estimate_tmle: family must be either "binomial" or "gaussian".')
  }

  # Because of quirk of program, delete observations with delta=0 if #>0
  # & < 10
  n = length(Y)
  inc = rep(T, n)

  # TODO: revisit this decision.
  if (!is.null(delta)) {
    num_missing = sum(delta == 0)
    if (num_missing > 0 && num_missing < 10) {
      inc[delta == 0] = F
    }
  }

  if (length(dim(W)) != 2) {
    stop("Error: W should have two dimensions. Instead its dimensions are:", paste(dim(W)), "\n")
  }

  Y = Y[inc]
  A = A[inc]
  W = W[inc, , drop = F]

  delta = delta[inc]

  # Check for any remaining missing data.
  # It is technically ok for Y or A to include missingness, because the
  # delta estimation is intended to estimate the missingness mechanism.
  # This needs more testing though.

  missing_vals = sum(is.na(W))
  if (missing_vals != 0) {
    cat("Warning: found", missing_vals, "NAs in W.\n")

    na_sum = sapply(W, function(col) sum(is.na(col)))
    cat("Columns with NAs:\n")
    print(na_sum[na_sum > 0])
  }

  missing_vals = sum(is.na(A))
  if (missing_vals != 0) {
    cat("Warning: found", missing_vals, "NAs in A.\n")
  }

  missing_vals = sum(is.na(Y[delta]))
  if (missing_vals != 0) {
    cat("Warning: found", missing_vals, "NAs in Y.\n")
  }

  # Here we are using tmle but not using the treatment effect estimate.
  # We're actually using the underlying variables to estimate Y_a.
  # TODO: disable this call, we're just running it now to double-check
  # the custom results.
  if (F) {
    tmle.1 = tmle::tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                      Q.SL.library = Q.lib, family = family, verbose = verbose)
  } else {
    tmle.1 = NULL
  }

  if (verbose) {
    cat("Estimating g. A distribution:\n")
    print(table(A))
  }

  min_g_cell = min(table(A))
  g_V = V
  if (min_g_cell < V) {
    g_V = max(min_g_cell, 2)
    if (verbose) {
      cat("A's minimum cell sizes", min_g_cell, "is less than the number of CV",
          "folds", V, "\n. Reducing folds to", g_V,  "\n")
    }
  }

  # Using modified version of tmle::estimateG
  #g_model = SuperLearner::SuperLearner(Y = A, X = W, SL.library = g.lib,
  #g = tmle_estimate_g(d = cbind(A[delta == 1], W[delta == 1, ]),
  g = tmle_estimate_g(d = cbind(A, W),
                      SL.library = g.lib,
                      verbose = verbose,
                      V = g_V,
                      outcome = "A",
                      message = "g")

  # Handle gBounds - code from tmle::tmle().
  if (length(gbound) == 1) {
    if (length(unique(A)) == 1) {
      # EY1 only, no controlled direct effect
      gbound = c(gbound, 1)
    } else {
      gbound = c(gbound, 1 - gbound)
    }
  }
  g$bound = gbound

  # Propensity score for treatment.
  #g1 = tmle.1$g$g1W
  g1 = g$g1W

  # Check if these two are highly correlated.
  if (!is.null(tmle.1) && verbose) {
    # This will generate a warning of SD is zero for either vector.
    # If this is the case we'll see an NA here.
    suppressWarnings(cat("Correlation of custom g to tmle-based g:",
        stats::cor(tmle.1$g$g1W, g1), "\n"))
  }

  # This is copied from within tmle::tmle()
  map_to_ystar = fluctuation == "logistic"

  if (verbose) cat("TMLE init stage1\n")

  # Run tmle stage 1 - this is basically just bounding & transforming Y.
  stage1 = tmle_init_stage1(Y = Y, Q = NULL,
                             A = A,
                             Delta = delta,
                             alpha = alpha,
                             Qbounds = Qbounds,
                             maptoYstar = map_to_ystar,
                             family = family)

  if (verbose) cat("TMLE q\n")

  # Estimate Qinit
  # cvQinit = F by default, meaning that we don't need CV.SuperLearner.
  q = tmle_estimate_q(Y = stage1$Ystar,
                      A = A,
                      W = W,
                      Q = stage1$Q,
                      Delta = delta,
                      SL.library = Q.lib,
                      family = family,
                      V = V,
                      verbose = verbose,
                      maptoYstar = map_to_ystar,
                      Qbounds = stage1$Qbounds)

  # Convert from a matrix to a df, so we can use $ to access elements.
  # Nevermind, this breaks the plogis code.
  # q$Q = data.frame(q$Q)

  # Check for all NaN in QAW.
  stopifnot(mean(is.nan(q$Q[, "QAW"])) == 0)

  # cat(class(q$Q), paste(dim(q$Q)), paste(colnames(q$Q)), "\n")
  # TODO: check if our custom q$QAW equals the tmle Q.

  # Specify random arguments from tmle::tmle()
  pDelta1 = NULL
  g.Deltaform = NULL

  # From tmle::tmle()
  ############################################
  if (verbose) cat("Estimating g.Delta (missingness mechanism)\n")

  g.z <- NULL
  g.z$type="No intermediate variable"
  g.z$coef=NA
  g.Delta <- suppressWarnings({
    tmle_estimate_g(d = data.frame(delta, Z=1, A, W),
                    pDelta1,
                    g.Deltaform,
                    g.lib,
                    id = id, V = V,
                    verbose = verbose,
                    message = "missingness mechanism",
                    outcome="D")
  })
  g1W.total <- .bound(g$g1W*g.Delta$g1W[,"Z0A1"], gbound)
  if (sum(is.na(g1W.total)) > 0) {
    if (verbose) {
      cat("Error, g1W.total has NAs:", sum(is.na(g1W.total)), "\n")
    }
  }
  g0W.total <- .bound((1-g$g1W)*g.Delta$g1W[,"Z0A0"], gbound)
  if(all(g1W.total==0)){g1W.total <- rep(10^-9, length(g1W.total))}
  if(all(g0W.total==0)){g0W.total <- rep(10^-9, length(g0W.total))}
  H1W <- A/g1W.total
  H0W <- (1-A)/g0W.total

  if (verbose) cat("Estimating epsilon\n")

  suppressWarnings(
    epsilon <- coef(glm(stage1$Ystar~-1 + offset(q$Q[,"QAW"]) + H0W + H1W, family=q$family, subset=delta==1))
  )
  epsilon[is.na(epsilon)] <- 0  # needed for EY1 calculation
  Qstar <- q$Q + c((epsilon[1]*H0W + epsilon[2]*H1W), epsilon[1]/g0W.total, epsilon[2]/g1W.total)
  colnames(Qstar) <- c("QAW", "Q0W", "Q1W")
  Ystar <- stage1$Ystar
  if (map_to_ystar) {
    Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1]
    Qstar <- plogis(Qstar)*diff(stage1$ab)+stage1$ab[1]
    q$Q <- plogis(q$Q)*diff(stage1$ab)+stage1$ab[1]
    Ystar <- Ystar*diff(stage1$ab)+stage1$ab[1]
  }
  colnames(q$Q) <- c("QAW", "Q0W", "Q1W")
  q$Q <- q$Q[,-1]

  if (verbose) cat("tmle::calcParameters\n")
  res <- tmle::calcParameters(Ystar, A, I.Z=rep(1, length(Ystar)), delta, g1W.total, g0W.total, Qstar,
                        mu1=mean(Qstar[,"Q1W"]), mu0=mean(Qstar[,"Q0W"]), id, family, 
                        obsWeights=rep(1, length(Ystar)))

  #returnVal <- list(estimates=res, Qinit=Q, g=g, g.Z=g.z, g.Delta=g.Delta, Qstar=Qstar[,-1], epsilon=epsilon)
  #class(returnVal) <- "tmle"
  ############################################




  # Calculate Qstar
  #QbarAW_star = plogis(qlogis(q$Q$QAW) + epsilon * h_aw)
  #Qbar1W_star = plogis(qlogis(q$Q$Q1W) + epsilon * h_1w)
  #Qbar0W_star = plogis(qlogis(q$Q$Q0W) + epsilon * h_0w)

  # Unit's estimated outcome under treatment: hat(Y) | A = 1, W
  # This is after the fluctuation step, so it is targeted.
  # Qst = tmle.1$Qstar[, 2]
  # Qst = tmle.1$Qstar$Q1W
  # Qst = Qbar1W_star
  Qst = Qstar[, "Q1W"]

  # E[Y | A = 1, W]
  theta = mean(Qst)

  # Influence curve
  IC = (A / g1) * (Y - Qst) + Qst - theta

  # Compile results.
  result = list(theta = theta, IC = IC, g_model = g$model, q_model = q$model,
                tmle = tmle.1, alpha = alpha,
                Qbounds = Qbounds,
                stage1_Qbounds = stage1$Qbounds,
                gbounds = g$bound,
                V = V,
                # May have had to reduce the # of SL folds for g, due to sparsity
                # in A.
                g_V = g_V,
                fluctuation = fluctuation,
                map_to_ystar = map_to_ystar, ab = stage1$ab)

  return(result)
}

