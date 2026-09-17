# CK: this is copied from TMLE package, with some modifications.
# Edits:
# - full superLearner model is returned.
# - set default for id argument.
# - dropped the branches varimpact cannot reach: the g1W and gform escape
#   hatches for externally supplied values, the outcome == "Z" intermediate
#   variable, and the fallback for SuperLearner < 2.0.
#
#-----------estimateG----------------
# Estimate factors of g
# 		P(A=1|W), P(Delta=1|Z,A,W)
# d - dataframe (A,W) or (Delta,Z,A,W)
# SL.library - algorithms to use for super learner estimation
#   id - subject identifier
# verbose - flag, whether or not to print messages
# message - printed when verbose=TRUE
# outcome - "A" for treatment, "D" for Delta (missingness)
# min_cell_size - smallest number of observations the rarer outcome class may
#           have before the covariate-adjusted fit is abandoned in favor of the
#           marginal proportion. 0 (the default) never abandons it.
# newdata - optional values to predict on (needed by tmleMSM function)
# d = [A,W] for treatment
# d = [Delta, Z,A,W] for missingness
#----------------------------------------
tmle_estimate_g <-
  function (d,
            SL.library,
            id=1:nrow(d),
            V = 10,
            stratify = TRUE,
            # TODO: consider changing default to method.NNloglik
            method = "method.NNLS",
            verbose = F,
            message = "",
            outcome = "A",
            min_cell_size = 0,
            newdata=d)  {

  cvControl = SuperLearner::SuperLearner.CV.control(V = V,
                                                    stratifyCV = stratify,
                                                    shuffle = TRUE)
  SL.ok <- FALSE
  m <- NULL
  coef <- NA
  type <- NULL
  if(verbose){cat("\tEstimating", message, "\n")}
  if (length(unique(d[,1]))==1) {
    g1W <- rep(1,nrow(d))
    type <- paste("No", strsplit(message, " ")[[1]][1])
    if (outcome == "D") {
      g1W <- cbind(Z0A0=g1W, Z0A1=g1W, Z1A0=g1W, Z1A1=g1W)
    }
  } else if (min(table(d[, 1])) < min_cell_size) {
    # The rarer class is too small to support a covariate-adjusted fit: a
    # V-fold cross-validation needs at least a couple of its observations in
    # every training split, and the learners resample again on top of that.
    # Below that floor SuperLearner does not fail cleanly. Each learner that
    # cannot fit the split errors, SuperLearner prints the error and gives the
    # learner weight 0, and if every learner fails the fallback is a glm that
    # separates perfectly. The ensemble that survives is the intercept, so fit
    # the intercept directly and say so.
    p <- mean(d[, 1])
    m <- marginal_fit(p)
    type <- "marginal"
    if (verbose) {
      cat("\tRarest class of", colnames(d)[1], "has", min(table(d[, 1])),
          "observations, fewer than", min_cell_size,
          "\n\tUsing the marginal proportion", signif(p, 4),
          "instead of a covariate-adjusted fit\n")
    }
    g1W <- rep(p, nrow(newdata))
    if (outcome == "D") {
      g1W <- cbind(Z0A0 = g1W, Z0A1 = g1W, Z1A0 = g1W, Z1A1 = g1W)
    }
  } else {
    SL.ok <- TRUE

    # With no adjustment variables - which happens when the analyzed
    # variable is the only column in the data - the conditional probability
    # reduces to the marginal one. Learners that regress on the covariates
    # cannot fit at all: SL.glm's "Y ~ ." errors on a zero-column data frame
    # ("'.' in formula and no 'data' argument"), so SuperLearner drops it
    # and falls back to SL.mean, which is that same marginal. Ask for it
    # directly, so the estimate is unchanged but the failed fits are not.
    g_library <- SL.library
    if (ncol(d) <= 1L) {
      if (verbose) {
        cat("\tNo adjustment variables; estimating", message,
            "with SL.mean alone.\n")
      }
      g_library <- "SL.mean"
    }

    sl_id = id
    # If IDs are unique we don't need to pass them to SL, which breaks
    # stratification.
    if (length(unique(id)) == length(id)) {
      sl_id = NULL
    }
    arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE],
                    family="binomial", method = method, SL.library=g_library, cvControl = cvControl, id = sl_id)
    # TODO: are we sure we want to suppress warnings here?
    suppressWarnings({
      # Suppress any package startup messages if we can.
      suppressPackageStartupMessages({
        m <- try(do.call(SuperLearner::SuperLearner, arglist))
      })
      # Set call to null because do.call() messes up that element.
      m$call = NULL
    })
    if(identical(class(m),"SuperLearner")) {
      g1W <- as.vector(m$SL.predict)
    } else {
      SL.ok <- FALSE
      cat("Error estimating g using SuperLearner. Defaulting to glm\n")
    }
    if (!SL.ok){
      if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
      form <- paste(paste(colnames(d)[1],"~1"), paste(colnames(d)[-1], collapse = "+"), sep="+")
      m <- glm(form, data=d, family="binomial")
      g1W <- predict(m, newdata=newdata, type="response")
      coef <- coef(m)
    }
    # Get counterfactual predicted values
    if (outcome == "D") {
      if(identical(class(m),"SuperLearner")){
        g1W <- cbind(predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response",
                             X=d[,-1,drop=FALSE], Y=d[,1])[[1]],
                     predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response",
                             X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
                     predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response",
                             X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
                     predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response",
                             X=d[,-1, drop=FALSE], Y=d[,1])[[1]])
      } else{
        g1W <- cbind(predict(m, newdata=data.frame(Z=0, A=0, newdata[,-(1:3), drop=FALSE]), type="response"),
                     predict(m, newdata=data.frame(Z=0, A=1, newdata[,-(1:3), drop=FALSE]), type="response"),
                     predict(m, newdata=data.frame(Z=1, A=0, newdata[,-(1:3), drop=FALSE]), type="response"),	     	     		  predict(m, newdata=data.frame(Z=1, A=1, newdata[,-(1:3), drop=FALSE]), type="response"))
      }
      colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")
    }
  }
  if(is.null(type)){ type <- class(m)[1]}
  returnVal <- list(g1W=g1W, coef=coef, type=type, model = m)
  if(type=="SuperLearner"){
    returnVal$SL.library=SL.library
    returnVal$coef=m$coef
  }
  return(returnVal)
}
