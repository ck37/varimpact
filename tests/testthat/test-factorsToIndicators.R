library(varimpact)
library(testthat)

context("factorsToIndicators")

# Create test dataset.

set.seed(1, "L'Ecuyer-CMRG")
N = 200

num_normal = 7
X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))

# Add some missing data to X.
miss_num = 10
for (i in 1:miss_num) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

X_fac = data.frame(lapply(1:ncol(X), FUN=function(col_i) as.factor(floor(abs(pmin(pmax(X[, col_i], -1), 1)*3)))))
dim(X_fac)
colnames(X_fac) = paste0("fac_", 1:ncol(X_fac))
colnames(X_fac)
summary(X_fac)

#########################

table(X_fac[, 1], useNA="ifany")

# Test a single factor.
# Use column 3 which has 3 missing values based on the seed
results = factors_to_indicators(X_fac[, 3, drop = F], verbose = T)
dim(results$data)
# We should have indicators for 1, 2, 3.
colnames(results$data)

dim(results$missing_indicators)
colnames(results$missing_indicators)

# We should have 1 missing data indicator, with 3 0s and the rest 1s.
table(results$missing_indicators[, 1])
expect_equal(min(table(results$missing_indicators[, 1])), 3)

expect_gt(ncol(results$missing_indicators), 0)

# Test multiple factors.
results = factors_to_indicators(X_fac[, 1:2, drop = F], verbose = T)

dim(results$data)
colnames(results$data)

dim(results$missing_indicators)
colnames(results$missing_indicators)

table(results$missing_indicators[, 1])

# We expect our missing indicator matrix to have more than 0 columns.
expect_is(results$missing_indicators, "matrix")

# Test the full factor dataframe.
results = factors_to_indicators(X_fac, verbose = T)
dim(results$data)
colnames(results$data)

dim(results$missing_indicators)
colnames(results$missing_indicators)
