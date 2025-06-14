# Make sure we're using the rebuilt package. Suppress any error if it isn't loaded.
try(detach(package:varimpact), silent = TRUE)
library(varimpact)
library(testthat)

# Create test dataset.
context("varimpact(). Dataset A: continuous variables")

# Set multicore-compatible seed.
set.seed(1, "L'Ecuyer-CMRG")
# Can't go below 90 without changing more varimpact default settings.
N = 100
num_normal = 5
X = data.frame(matrix(rnorm(N * num_normal), N, num_normal))

# Systematic Y generation.
Y = .2 * X[, 1] + 1 * X[, 2] - 0.8 * X[, 3] + .1 * X[, 3] * X[, 4] - .2 * abs(X[, 4])

# Binary distribution via the binomial.
Y_bin = rbinom(N, 1, plogis(Y))

# Gaussian distribution.
Y_gaus = Y + rnorm(N, 0, 1)

# Add some missing data to X.
miss_num = 10
for (i in 1:miss_num) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

# Basic test - binary outcome.
#future::plan("multiprocess")
future::plan("sequential")
vim = varimpact(Y = Y_bin, data = X[, 1:3], V = 3L,
                Q.library = c("SL.mean", "SL.glm"),
                g.library = c("SL.mean", "SL.glm"),
                verbose = TRUE,
                verbose_tmle = FALSE, bins_numeric = 3L)
# Takes 25 seconds.
vim$time
# Be explict about printing for code coverage of tests.
print(vim)
vim$results_all
vim$results_by_fold
# names(vim)
# exportLatex testing moved to test-exportLatex.R

# And try a gaussian outcome.
vim = varimpact(Y = Y_gaus, data = X[, 1:3], V = 3L, verbose = TRUE,
                family = "gaussian")
print(vim)

# Test imputation
vim = varimpact(Y = Y_bin, data = X[, 1:3], verbose = TRUE, impute = "zero")
vim = varimpact(Y = Y_bin, data = X[, 1:3], verbose = TRUE, impute = "median")
vim = varimpact(Y = Y_bin, data = X[, 1:4], verbose = TRUE, impute = "knn")

# Test a subset of columns using A_names.
vim = varimpact(Y = Y_bin, data = X, A_names = colnames(X)[1:2], verbose = TRUE)
print(vim)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Test parallelization
  future::plan("multiprocess", workers = 2)
  vim = varimpact(Y = Y_bin, data = X[, 1:3], verbose = TRUE)
  print(vim)
}

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Test parallelization via snow.
  cl = snow::makeCluster(2L)
  future::plan("cluster", workers = cl)
  vim = varimpact(Y = Y_bin, data = X[, 1:4], verbose = TRUE)
  vim
  snow::stopCluster(cl)
}

context("varimpact(). Dataset B: factor variables")

# Set a new multicore-compatible seed.
set.seed(2, "L'Ecuyer-CMRG")

X_fac = data.frame(lapply(1:ncol(X),
                          function(col_i)
                        as.factor(floor(abs(pmin(pmax(X[, col_i], -1), 1) * 3)))))
dim(X_fac)
colnames(X_fac) = paste0("fac_", 1:ncol(X_fac))
colnames(X_fac)
summary(X_fac)

# Return to sequential execution for now.
future::plan("sequential")

# Basic factor test.
vim = varimpact(Y = Y_bin, data = X_fac[, 1:3], V = 2L, verbose = TRUE)
print(vim)

# And gaussian
vim = varimpact(Y = Y_gaus, data = X_fac[, 1:3], V = 2L, verbose = TRUE,
                family = "gaussian")
print(vim)

# Only run in RStudio so that automated CRAN checks don't give errors.
# Disabled for now - need to review.
if (F && .Platform$GUI == "RStudio") {
  # Test parallelization.
  future::plan("multiprocess")

  # Try a snow cluster, which does return the output to STDOUT.
  if (F) {
    # Run manually when debugging.
    cores = RhpcBLASctl::get_num_cores()
    capture.output({ cl = snow::makeCluster(cores, type="SOCK", outfile = "") })
    doSNOW::registerDoSNOW(cl)
    parallel::setDefaultCluster(cl)
    foreach::getDoParName()
  }

  # Factor variables with parallelization.
  # TOFIX: This does not complete currently if fac_4 is included.
  # I think it is due to HOPACH never completing.
  vim = varimpact(Y = Y_bin, data = X_fac[, 1:3],
                  #A_names = c(_4", "fac_2"),
                  verbose = TRUE)
  vim

  # Return to single core usage.

  # Run manually when debugging, if the snow cluster was used.
  if (F) {
    ck37r::stop_cluster(cl)
  }

}

context("varimpact(). Dataset C: numeric and factor variables")

#################################
# Combined numeric and factor test.
X_combined = cbind(X[1:3], X_fac[4:5])

# Basic combined test.
vim = varimpact(Y = Y_bin, data = X_combined, V = 2, verbose = TRUE)
print(vim)

# And gaussian
vim = varimpact(Y = Y_gaus, data = X_combined, V = 2, verbose = TRUE,
                family = "gaussian")
print(vim)

context("varimpact() .Dataset D: basic example")

####################################
# Create test dataset.
set.seed(1, "L'Ecuyer-CMRG")

N = 100
num_normal = 7
X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y = rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
# Add some missing data to X so we can test imputation.
for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

####################################
# Basic example
# TODO: fix warnings here, due to failed folds.
# TOFIX: there is an error here on the numeric variables.
# task 3 failed - "attempt to select less than one element in get1index"
# X_3 seems to be causing the problem - need to investigate why.
vim = varimpact(Y = Y, data = X, A_names = colnames(X)[c(1, 2, 4:7)],
                verbose = TRUE, parallel = FALSE)
print(vim)
vim$results_all
vim$results_by_fold
# In this test all variables are significant, which is rare.
# exportLatex testing moved to test-exportLatex.R
