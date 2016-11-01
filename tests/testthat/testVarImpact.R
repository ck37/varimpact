# Make sure we're using the rebuilt package. Suppress any error if it isn't loaded.
tryCatch( detach(package:varImpact), error = function(e) invisible() )
library(varImpact)
library(testthat)

# Create test dataset.
context("varImpact(). Dataset A: continuous variables")

# Set multicore-compatible seed.
set.seed(1, "L'Ecuyer-CMRG")
N = 100
num_normal = 5
X = data.frame(matrix(rnorm(N * num_normal), N, num_normal))

# Systematic Y generation.
Y = .2 * X[, 1] + .1 * X[, 2] - .2 * X[, 3] + .1 * X[, 3] * X[, 4] - .2 * abs(X[, 4])

# Binary distribution via the binomial.
Y_bin = rbinom(N, 1, plogis(Y))

# Gaussian distribution.
Y_gaus = Y + rnorm(N, 0, 1)

# Add some missing data to X.
miss_num = 10
for (i in 1:miss_num) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

# Basic test - binary outcome.
vim = varImpact(Y = Y_bin, data = X[, 1:3], V = 2, verbose=T)
vim$time
# Be explict about printing for code coverage of tests.
print(vim)
vim$results_all
vim$results_by_fold
# names(vim)
exportLatex(vim)
# Clean up - will get a warning if there were no consistent results.
suppressWarnings({
  file.remove(c("varimpByV.tex", "varImpAll.tex", "varimpConsistent.tex"))
})

# And try a gaussian outcome.
vim = varImpact(Y = Y_gaus, data = X[, 1:4], V = 2, verbose=T, family="gaussian")
vim

# Test imputation
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose = T, impute = "zero")
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose = T, impute = "median")
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose = T, impute = "knn")


# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Test parallelization via doMC.
  doMC::registerDoMC()
  # Check how many cores we're using.
  foreach::getDoParWorkers()
  vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose = T)
}

# Test disabling parallelization.
vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose = T, parallel = F)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Return to single core usage.
  foreach::registerDoSEQ()
  foreach::getDoParWorkers()
}

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Test parallelization via doSnow.
  cluster = snow::makeCluster(2)
  doSNOW::registerDoSNOW(cluster)
  # Check that we're using the snow cluster.
  expect_equal(foreach::getDoParName(), "doSNOW")
  expect_equal(foreach::getDoParWorkers(), 2)

  vim = varImpact(Y = Y_bin, data = X[, 1:4], verbose=T)
  vim
}

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  snow::stopCluster(cluster)
  # Return to single core usage.
  foreach::registerDoSEQ()
}

context("varImpact(). Dataset B: factor variables")

# Set a new multicore-compatible seed.
set.seed(2, "L'Ecuyer-CMRG")

X_fac = data.frame(lapply(1:ncol(X), FUN=function(col_i)
  as.factor(floor(abs(pmin(pmax(X[, col_i], -1), 1)*3)))))
dim(X_fac)
colnames(X_fac) = paste0("fac_", 1:ncol(X_fac))
colnames(X_fac)
summary(X_fac)

# Basic factor test.
# TODO: this generates multiple errors for fac_4
vim = varImpact(Y = Y_bin, data = X_fac[, 1:3], V = 2, verbose=T)
vim
# And gaussian
vim = varImpact(Y = Y_gaus, data = X_fac[, 1:3], V = 2, verbose=T, family="gaussian")
vim

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Test parallelization.
  doMC::registerDoMC()
  # Check how many cores we're using.
  foreach::getDoParWorkers()

  # Factor variables with parallelization.
  vim = varImpact(Y = Y_bin, data = X_fac[, 1:4], verbose=T)
  vim
}

context("varImpact(). Dataset C: numeric and factor variables")

#################################
# Combined numeric and factor test.
X_combined = cbind(X[1:3], X_fac[4:5])

# Basic combined test.
vim = varImpact(Y = Y_bin, data = X_combined, V = 2, verbose=T)
vim

# And gaussian
vim = varImpact(Y = Y_gaus, data = X_combined, V = 2, verbose=T, family="gaussian")
vim

#################################
# mlbench BreastCancer dataset.

context("BreastCancer dataset")

data(BreastCancer, package="mlbench")
data = BreastCancer

set.seed(3, "L'Ecuyer-CMRG")

# Reduce to a dataset of 100 observations to speed up testing.
data = data[sample(nrow(data), 100), ]

# Create a numeric outcome variable.
data$Y = as.numeric(data$Class == "malignant")
table(data$Y)

X = subset(data, select=-c(Y, Class, Id))
dim(X)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  doMC::registerDoMC()
}
# This takes 1-2 minutes.
vim = varImpact(Y = data$Y, X, verbose = T)
vim$time
vim

context("varImpact() .Dataset D: basic example")

####################################
# Create test dataset.
set.seed(1, "L'Ecuyer-CMRG")

N <- 100
num_normal <- 7
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
# Add some missing data to X so we can test imputation.
for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA

####################################
# Basic example
# TODO: fix warnings here, due to failed folds.
vim <- varImpact(Y = Y, data = X, verbose = T)
vim
vim$results_all
vim$results_by_fold
# In this test all variables are significant, which is rare.
exportLatex(vim)
# Clean up
suppressWarnings({ # Suppress a warning when no results are consistent.
  file.remove(c("varimpByV.tex", "varImpAll.tex", "varimpConsistent.tex"))
})
