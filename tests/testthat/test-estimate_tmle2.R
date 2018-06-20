library(varimpact)
library(SuperLearner)
library(tmle)
library(testthat)

# Create test dataset.
context("Dataset A: continuous variables")

# Set multicore-compatible seed.
set.seed(1, "L'Ecuyer-CMRG")

# Simulation code from PH 295, Fall 2016 --
# https://github.com/wilsoncai1992/PH295-lab/blob/master/lab3/Lab3_Lecture.Rmd

# structural equation for W_1
# takes as input a vector U_W1 and returns a vector evaluating
# f_{W,1}(U_W1)
f_W1 <- function(U_W1){
  return(U_W1)
}

# structural equation for W_2
# takes as input a vector U_W2 and returns a vector evaluating
# f_{W,2}(U_W2)
f_W2 <- function(U_W2){
  return(U_W2)
}

# structural equation for A
f_A <- function(W_1, W_2, U_A){
  return(as.numeric(plogis(W_1 - W_2 + U_A) > 0.5))
}

# structural equation for Y
f_Y <- function(W_1, W_2, A, U_Y){
  return(-W_1 + W_2 + A - U_Y)
}

# function to draw n observations from an scm
# n = the number of observations to draw
# returns a data.frame with named columns
simObsSCM <- function(n){
  ## first we draw the errors
  # draw U_{W,1}
  U_W1 <- rbinom(n,1,0.5)
  # draw U_{W,2}
  U_W2 <- rbinom(n,1,0.5)
  # draw U_A
  U_A <- rnorm(n,0,1)
  # draw U_Y
  U_Y <- rnorm(n,0,1)

  ## now we can evaluate the observations sequentially
  # evaluate W_1
  W_1 <- f_W1(U_W1)
  #evaluate W_2
  W_2 <- f_W2(U_W2)
  # evaluate A
  A <- f_A(W_1 = W_1, W_2 = W_2, U_A = U_A)
  # evaluate Y
  Y <- f_Y(W_1 = W_1, W_2 = W_2, A = A, U_Y = U_Y)

  ## return a data.frame object
  out <- data.frame(W_1 = W_1, W_2 = W_2, A = A, Y = Y)
  return(out)
}

data = simObsSCM(100)
summary(data)

sl_lib = c("SL.glmnet", "SL.glm", "SL.mean")

# Estimate g and Q
result = varimpact::estimate_tmle2(Y = data$Y, A = data$A,
  W = data[, c("W_1", "W_2")], family = "gaussian",
  Q.lib = sl_lib,
  g.lib = sl_lib,
  verbose = T)
