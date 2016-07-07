library(varImpact)
library(testthat)

# Create test dataset.
context("Dataset A: continuous variables")

set.seed(1)
N = 200

num_normal = 7
X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))

# Binary distribution via the binomial.
Y_bin = rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
# Gaussian distribution.
Y_gaus = .2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4]) + rnorm(N, 0, 1)

# Add some missing data to X.
miss_num = 10
for (i in 1:miss_num) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

# Basic test - binary outcome.
vim = varImpact(Y = Y_bin, data = X[, 1:4], V = 2, verbose=T)
vim
vim$results_all
# names(vim)
exportLatex(vim, dir = "results")

# And try a gaussian outcome.
vim = varImpact(Y = Y_gaus, data = X[, 1:4], V = 2, verbose=T, family="gaussian")

# Test imputation
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T, impute="zero")
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T, impute="median")
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T, impute="knn")


# Test parallelization via doMC.
doMC::registerDoMC()
# Check how many cores we're using.
foreach::getDoParWorkers()
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T)

# Test disabling parallelization.
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T, parallel = F)
# Return to single core usage.
foreach::registerDoSEQ()
foreach::getDoParWorkers()

# Test parallelization via doSnow.
cluster = snow::makeCluster(2)
doSNOW::registerDoSNOW(cluster)
# Check that we're using the snow cluster.
foreach::getDoParName()
foreach::getDoParWorkers()
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T)
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
vim = varImpact(Y = Y_bin, data = X_fac[, 1:4], V = 2, verbose=T)
# And gaussian
vim = varImpact(Y = Y_gaus, data = X_fac[, 1:4], V = 2, verbose=T, family="gaussian")

# Test parallelization.
doMC::registerDoMC()
# Check how many cores we're using.
foreach::getDoParWorkers()

# Factor variables with parallelization.
vim = varImpact(Y = Y_bin, data = X_fac[, 1:4], verbose=T)
vim

context("Dataset C: numeric and factor variables")

#################################
# Combined numeric and factor test.
X_combined = cbind(X[1:3], X_fac[5:7])

# Basic combined test.
vim = varImpact(Y = Y_bin, data = X_combined, V = 2, verbose=T)
vim
# And gaussian
vim = varImpact(Y = Y_gaus, data = X_combined, V = 2, verbose=T, family="gaussian")

#################################
# mlbench BreastCancer dataset.
data(BreastCancer, package="mlbench")
data = BreastCancer

# Create a numeric outcome variable.
data$Y = as.numeric(data$Class == "malignant")
table(data$Y)

# Use multicore parallelization to speed up processing.
doMC::registerDoMC()
# This takes 1-3 minutes.
vim = varImpact(Y = data$Y, data = subset(data, select=-c(Y, Class, Id)))
vim
