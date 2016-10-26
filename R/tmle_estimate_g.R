# CK: this is copied from TMLE package
# Edits:
# - full superLearner model is returned.
# - set default for id argument.
#
#-----------estimateG----------------
# Estimate factors of g
# 		P(A=1|W), P(Z=1|A,W), P(Delta=1|Z,A,W)
# d - dataframe (A,W), (Z,A,W), or (Delta,Z,A,W)
# g1W - optional vector/matrix  of externally estimated values
# gform - optionalformula to use for glm
# SL.library - algorithms to use for super learner estimation
#   id - subject identifier
# verbose - flag, whether or not to print messages
# message - printed when verbose=TRUE
# outcome - "A" for treatment, "Z" for intermediate variable,
#           "D" for Delta (missingness)
# newdata - optional values to predict on (needed by tmleMSM function)
# d = [A,W] for treatment
# d = [Z,A,W] for intermediate
# d = [Delta, Z,A,W for missingness]
#----------------------------------------
#' @importFrom stats packageDescription
tmle_estimate_g <- function (d,g1W = NULL, gform = NULL,SL.library, id=1:nrow(d), verbose = F, message = "", outcome="A", newdata=d)  {
  SL.version <- 2
  SL.ok <- FALSE
  m <- NULL
  coef <- NA
  type <- NULL
  if (is.null(g1W)){
    if(verbose){cat("\tEstimating", message, "\n")}
    if (length(unique(d[,1]))==1) {
      g1W <- rep(1,nrow(d))
      type <- paste("No", strsplit(message, " ")[[1]][1])
      if(outcome=="Z"){
        g1W <- cbind(A0=g1W, A1=g1W)
      } else if (outcome=="D"){
        g1W <- cbind(Z0A0=g1W, Z0A1=g1W, Z1A0=g1W, Z1A1=g1W)
      }
    } else {
      if (is.null(gform)){
        SL.ok <- TRUE
        old.SL <- packageDescription("SuperLearner")$Version < SL.version
        if(old.SL){
          arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE], family="binomial", SL.library=SL.library, V=5, id=id)
        } else {
          arglist <- list(Y=d[,1], X=d[,-1, drop=FALSE], newX=newdata[,-1, drop=FALSE], family="binomial", SL.library=SL.library, cvControl=list(V=5), id=id)
        }
        suppressWarnings({
          m <- try(do.call(SuperLearner::SuperLearner, arglist))
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
      } else {
        form <- try(as.formula(gform))
        if(class(form)== "formula") {
          m <- try(glm(form,  data=d, family="binomial"))
          if (class(m)[1]=="try-error"){
            if(verbose){cat("\tInvalid formula supplied. Running glm using main terms\n")}
            form <- paste(colnames(d)[1],"~1 + ", paste(colnames(d)[-1], collapse = "+"), sep="")
            m <- glm(form, data=d, family="binomial")
          } else {
            type <- "user-supplied regression formula"
          }
        } else {
          if(verbose){cat("\tRunning main terms regression for 'g' using glm\n")}
          form <- paste(colnames(d)[1],"~1", paste(colnames(d)[-1], collapse = "+"), sep="+")
          m <- glm(form, data=d, family="binomial")
        }
        g1W <- predict(m, newdata=newdata, type="response")
        coef <- coef(m)
      }
      # Get counterfactual predicted values
      if(outcome=="Z"){
        if(identical(class(m),"SuperLearner")){
          g1W <- cbind(predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response",
                               X=d[,-1, drop=FALSE], Y=d[,1])[[1]],
                       predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response",
                               X=d[,-1, drop=FALSE], Y=newdata[,1])[[1]])
        } else {
          g1W <- cbind(predict(m, newdata=data.frame(A=0, newdata[,-(1:2), drop=FALSE]), type="response"),
                       predict(m, newdata=data.frame(A=1, newdata[,-(1:2), drop=FALSE]), type="response"))
        }
        colnames(g1W) <- c("A0", "A1")

      } else if (outcome=="D"){
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
  } else {
    type <- "user-supplied values"
    if(outcome=="Z") {
      colnames(g1W) <- c("A0", "A1")
    } else if (outcome=="D"){
      colnames(g1W) <- c("Z0A0", "Z0A1", "Z1A0", "Z1A1")[1:ncol(g1W)]
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
