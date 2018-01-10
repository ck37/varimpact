library(testthat)
library(varImpact)

#################################
# mlbench BreastCancer dataset.

context("BreastCancer dataset")

data(BreastCancer, package = "mlbench")
data = BreastCancer

set.seed(3, "L'Ecuyer-CMRG")

# Reduce to a dataset of 200 observations to speed up testing.
data = data[sample(nrow(data), 200), ]

# Create a numeric outcome variable.
data$Y = as.numeric(data$Class == "malignant")
table(data$Y)

X = subset(data, select = -c(Y, Class, Id))
dim(X)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  future::plan("multiprocess", workers = 2)
}
# This takes 1-2 minutes.
vim = varImpact(Y = data$Y, X, verbose = T, verbose_tmle = F)
vim$time
vim

# Test a subset of columns for A_names.
colnames(X)[1:3]
vim = varImpact(Y = data$Y, X, A_names = colnames(X)[1:3], verbose = T)
vim$time
vim
vim$results_all

# Return to single core usage.
future::plan("sequential")
