
library(testthat)

# Create test dataset.
context("Dataset A: continuous variables")

set.seed(1)
N = 200

num_normal = 7
X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))

Y = rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))

# Add some missing data to X.
miss_num = 10
for (i in 1:5) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

# This works.
vim = varImpact(Y = Y, data = X, V = 2, verbose=T)
vim
vim$results_all
# names(vim)
exportLatex(vim, dir = "results")


# Test imputation
vim = varImpact(Y = Y, data = X, V = 2, verbose=T, impute="median")
vim = varImpact(Y = Y, data = X, V = 2, verbose=T, impute="knn")

#context("Dataset B: factor variables")
#context("Dataset C: numeric and factor variables")
