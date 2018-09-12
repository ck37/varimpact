library(testthat)
library(varimpact)

#################################
# mlbench BreastCancer dataset.

context("BreastCancer dataset")

data(BreastCancer, package = "mlbench")
data = BreastCancer
names(data) = tolower(names(data))

set.seed(3, "L'Ecuyer-CMRG")

# Reduce to a dataset of 200 observations to speed up testing.
data = data[sample(nrow(data), 200), ]

# Create a numeric outcome variable.
data$y = as.numeric(data$class == "malignant")
table(data$y)

x = subset(data, select = -c(y, class, id))
dim(x)

# Only run in RStudio so that automated CRAN checks don't give errors.
if (.Platform$GUI == "RStudio") {
  # Use multicore parallelization to speed up processing.
  future::plan("multiprocess", workers = 2)
}
# This takes 1-2 minutes.
vim = varimpact(Y = data$y, x, verbose = TRUE, verbose_tmle = FALSE)
vim$time
vim

# Test a subset of columns for A_names.
colnames(x)[1:3]
vim = varimpact(Y = data$y, x, A_names = colnames(x)[1:3], verbose = TRUE)
vim$time
vim
vim$results_all

# Return to single core usage.
future::plan("sequential")
