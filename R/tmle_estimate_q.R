# CK:
# Edits:
#  - SuperLearner needs to have saveFitLibrary = T. (critical tweak)
#  - return Q model.
#  - Stop if Qbounds is null.
#  - add default for id argument.
#  - add NULL default for Z, Q, QForm arguments.
#-----------estimateQ----------------
#' purpose: estimate Q=E(Y |Z, A,W) data-adaptively,
#'
#' unless super learner not available, or user specifies
#' initial values or a regression formula
#' arguments:
#' @param Y outcome
#' @param A treatment indicator (1=treatment, 0=control)
#' @param W baseline covariates
#' @param	Delta missingness indicator
#' @param	Q optional externally estimated values for Q
#' @param	Qbounds bounds for predicted values
#' @param	Qform optional regression formula to use for glm if
#	        non-data adaptive estimation specified
#' @param maptoYstar if TRUE, using logistic fluctuation for bounded, continuous outcomes
# 		estimation inital Q on linear scale, bounded by (0,1),and return on logit scale
#		(will work if family=poisson)
#' @param	SL.library library of prediction algorithms for Super Learner
#' @param family regression family
#' @param	id subject identifier
#' @param V number of folds for SuperLearner
#' @param verbose Set T for extra output
# returns matrix of linear predictors for Q(A,W), Q(0,W), Q(1,W),
#		family for stage 2 targeting
#		coef, NA, unless Q is estimated using a parametric model
# 		type, estimation method for Q
#----------------------------------------
#' @export
tmle_estimate_q <-
  function(Y,
           A,
           W,
           Delta,
           Q = NULL,
           Qbounds,
           Qform = NULL,
           maptoYstar,
           SL.library,
           family,
           id = 1:length(Y),
           V = 10,
           verbose = F) {

  if (is.null(Qbounds)) stop("Qbounds must be defined.")
  Qfamily <- family
  m <- NULL
  coef <- NA
  type <- "user-supplied values"
  if(is.null(Q)){
    if(verbose) { cat("\tEstimating initial regression of Y on A and W\n")}
    Q <- matrix(NA, nrow=length(Y), ncol = 3)
    colnames(Q)<- c("QAW", "Q0W", "Q1W")
    if(!(is.null(Qform))){
      if(identical(as.character(as.formula(Qform)), c("~","Y", "."))){
        Qform <- paste("Y~A+", paste(colnames(W), collapse="+"))
      }
      m <- suppressWarnings(glm(Qform, data=data.frame(Y,A,W, Delta), family=family, subset=Delta==1))
      Q[,"QAW"] <- predict(m, newdata=data.frame(Y,A,W), type="response")
      Q[,"Q0W"] <- predict(m, newdata=data.frame(Y,A=0,W), type="response")
      Q[,"Q1W"] <- predict(m, newdata=data.frame(Y,A=1,W), type="response")
      coef <- coef(m)
      type="glm, user-supplied model"
    } else {
      if(verbose) {cat("\t using SuperLearner\n")}
      n <- length(Y)
      X <- data.frame(A,W)
      X00 <- data.frame(A=0, W)
      X01 <- data.frame(A=1, W)
      newX <- rbind(X, X00, X01)
      arglist <- list(Y=Y[Delta==1],X=X[Delta==1, , drop = FALSE], newX=newX, SL.library=SL.library,
                      cvControl=list(V=V), family=family, control = list(saveFitLibrary=T), id=id[Delta==1])
      suppressWarnings({
        # CK: try to eliminate messages from loading packages.
        out = utils::capture.output({
          suppressPackageStartupMessages({
            m <- try(do.call(SuperLearner::SuperLearner, arglist))
          })
        })
        # Set call to null because do.call() messes up that element.
        m$call = NULL
      })
      if (identical(class(m),"SuperLearner")){
        #if (verbose) print(m)
        Q[,"QAW"] <- m$SL.predict[1:n]
        Q[,"Q0W"] <- m$SL.predict[(n+1):(2*n)]
        Q[,"Q1W"] <- m$SL.predict[(2*n+1):(3*n)]
        type <- "SuperLearner"
      } else {
        stop("Super Learner failed when estimating Q. Exiting program\n")
      }
    }
  }
  if(is.na(Q[1,1]) | identical(class(m), "try-error")){
    if(verbose) {cat("\t Running main terms regression for 'Q' using glm\n")}
    Qform <- paste("Y~A+", paste(colnames(W), collapse="+"))
    m <- glm(Qform, data=data.frame(Y,A,W, Delta), family=family, subset=Delta==1)
    Q[,"QAW"] <- predict(m, newdata=data.frame(Y,A,W), type="response")
    Q[,"Q1W"] <- predict(m, newdata=data.frame(Y,A=1,W), type="response")
    Q[,"Q0W"] <- predict(m, newdata=data.frame(Y,A=0,W), type="response")
    coef <- coef(m)
    type="glm, main terms model"
  }
  Q <- varimpact::.bound(Q, Qbounds)
  if(maptoYstar | identical(Qfamily,"binomial") | identical(Qfamily, binomial)){
    Q <- qlogis(Q)
    Qfamily <- "binomial"
  } else if (identical(Qfamily, "poisson") | identical(Qfamily, poisson)) {
    Q <- log(Q)
    Qfamily <- "poisson"
  }
  Qinit <- list(Q=Q, family=Qfamily, coef=coef, type=type, model = m)
  if(type=="SuperLearner"){
    Qinit$SL.library=SL.library
    Qinit$coef=m$coef
  }
  return(Qinit)
}
