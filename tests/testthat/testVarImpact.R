
library(testthat)

# Create test dataset.
context("Dataset A - continuous only variables")

set.seed(1)
N = 2000

num_normal = 10
X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))

Y = rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))

# This works.
vim = varImpact(Y = Y, data = X, V = 2, verbose=T)
vim
vim$results_all
names(vim)
exportLatex(vim, dir = "results")
