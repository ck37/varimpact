# Copied from tmle package.
#---------- function .initStage1 ---------------
# Bound Y, map to Ystar if applicable, and
# set boundson on Q and enforce on user-specified values
# returns
#   Ystar - outcome values (between [0,1] if maptoYstar=TRUE)
#   Q - matrix of user-specified values
#   Qbounds - bounds on predicted values for Q (10% wider at each end then
# 			observed range of Y
#			(-Inf,+Inf) is default for linear regression
#   ab - bounding levels used to transform Y to Ystar
#-----------------------------------------------
tmle_init_stage1 <- function(Y,A, Q, Q.Z1=NULL, Delta, Qbounds, alpha, maptoYstar, family){
  if(family=="binomial") {Qbounds <- c(0,1)}
  if(is.null(Qbounds)) {
    if(maptoYstar){
      Qbounds <- range(Y[Delta==1])
      Qbounds <- Qbounds + .1*c(-abs(Qbounds[1]),abs(Qbounds[2]))
    } else {
      Qbounds <- c(-Inf, Inf)
    }
  }
  if(!is.null(Q)){
    QAW <- (1-A)*Q[,1] + A*Q[,2]
    Q <- cbind(QAW, Q0W=Q[,1], Q1W=Q[,2])
  }
  if(!is.null(Q.Z1)){
    Q <- cbind(Q, Q0W.Z1=Q.Z1[,1], Q1W.Z1=Q.Z1[,2])
  }
  ab <- c(0,1)
  Ystar <- Y
  if(maptoYstar){
    Ystar <- .bound(Y, Qbounds)
    if(!is.null(Q)){
      Q <- .bound(Q, Qbounds)
    }
    if(0 >= alpha | 1 <= alpha){
      alpha <- .995
      warning(paste("\n\talpha must be between 0 and 1, alpha reset to",alpha,"\n"),
              immediate. = TRUE)
    }
    ab <- range(Ystar, na.rm=TRUE)
    Ystar[is.na(Ystar)] <- 0
    Ystar <- (Ystar-ab[1])/diff(ab)
    if(!is.null(Q)){Q <- (Q-ab[1])/diff(ab)}
    Qbounds <- c(alpha, 1-alpha)
  }
  return(list(Ystar=Ystar, Q=Q, Qbounds=Qbounds, ab=ab))
}
