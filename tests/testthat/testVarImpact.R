library(varImpact)
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

# Basic test.
vim = varImpact(Y = Y, data = X[, 1:4], V = 2, verbose=T)
vim
vim$results_all
# names(vim)
exportLatex(vim, dir = "results")


# Test imputation
vim = varImpact(Y = Y, data = X[, 1:4], verbose=T, impute="zero")
vim = varImpact(Y = Y, data = X[, 1:4], verbose=T, impute="median")
vim = varImpact(Y = Y, data = X[, 1:4], verbose=T, impute="knn")


# Test parallelization via doMC.
doMC::registerDoMC()
# Check how many cores we're using.
foreach::getDoParWorkers()
vim = varImpact(Y = Y, data = X[, 1:4], verbose=T)

# Test disabling parallelization.
vim = varImpact(Y = Y, data = X[, 1:4], verbose=T, parallel = F)
# Return to single core usage.
foreach::registerDoSEQ()
foreach::getDoParWorkers()

# Test parallelization via doSnow.
cluster = snow::makeCluster(2)
doSNOW::registerDoSNOW(cluster)
# Check that we're using the snow cluster.
foreach::getDoParName()
foreach::getDoParWorkers()
vim = varImpact(Y = Y, data = X[, 1:4], verbose=T)
vim
snow::stopCluster(cluster)
# Return to single core usage.
foreach::registerDoSEQ()

context("Dataset B: factor variables")

X_fac = data.frame(lapply(1:ncol(X), FUN=function(col_i) as.factor(floor(abs(pmin(pmax(X[, col_i], -1), 1)*3)))))
dim(X_fac)
colnames(X_fac) = paste0("fac_", 1:ncol(X_fac))
colnames(X_fac)
summary(X_fac)

# Basic factor test.
vim = varImpact(Y = Y, data = X_fac[, 1:4], V = 2, verbose=T)

# Test parallelization.
doMC::registerDoMC()
# Check how many cores we're using.
foreach::getDoParWorkers()

# Factor variables with parallelization.
vim = varImpact(Y = Y, data = X_fac[, 1:4], verbose=T)
vim

context("Dataset C: numeric and factor variables")

X_combined = cbind(X[1:3], X_fac[5:7])

# Basic factor test.
vim = varImpact(Y = Y, data = X_combined, V = 2, verbose=T)
vim
