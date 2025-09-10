# Minimal test to reproduce GitHub issue #32: Missing Y values bug
# This test directly tests the problematic code without requiring the full varimpact package

# Reproduce the exact simulation from the GitHub issue
set.seed(1, "L'Ecuyer-CMRG")
N <- 200
num_normal <- 4
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))

# Add some missing data to X so we can test imputation
for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA

# Add missing data to Y - exactly as in the GitHub issue
# This creates 11 missing Y values, which triggers the bug
Y[c(4,6,7,8,11,15,20,21,28,32,72)] <- NA

cat("Number of missing Y values:", sum(is.na(Y)), "\n")
cat("This should be 11, which is > 10 and triggers the bug\n")

# Simulate the problematic code from vim-factors.R lines 175-186
# This is what happens inside the varimpact function

# Create some dummy data to simulate the internal state
folds <- sample(1:2, N, replace = TRUE)  # 2-fold CV
fold_k <- 1

# Simulate training data (all data not in this fold)
Yt <- Y[folds != fold_k]
At <- as.factor(sample(c("A", "B", "C"), sum(folds != fold_k), replace = TRUE))  # Dummy factor variable
Wtsht <- matrix(rnorm(sum(folds != fold_k) * 3), ncol = 3)  # Dummy adjustment variables

# This is the key line from vim-factors.R line 175
deltat <- as.numeric(!is.na(Yt) & !is.na(At))

cat("Number of missing observations (deltat == 0):", sum(deltat == 0), "\n")

# This is the problematic code from lines 181-186
cat("Testing the problematic condition: sum(deltat == 0) < 10\n")
cat("sum(deltat == 0) =", sum(deltat == 0), "\n")
cat("sum(deltat == 0) < 10 =", sum(deltat == 0) < 10, "\n")

if (sum(deltat == 0) < 10) {
  cat("CLEANUP CODE RUNS: Removing missing observations\n")
  Yt_original_length <- length(Yt)
  Yt <- Yt[deltat == 1]
  At <- At[deltat == 1]
  Wtsht <- Wtsht[deltat == 1, , drop = FALSE]
  deltat <- deltat[deltat == 1]
  cat("Yt length before cleanup:", Yt_original_length, "after cleanup:", length(Yt), "\n")
} else {
  cat("CLEANUP CODE DOES NOT RUN: Missing observations remain in data\n")
  cat("This is the bug! TMLE will receive data with missing values and fail\n")
}

cat("Final Yt length:", length(Yt), "\n")
cat("Final number of missing Yt values:", sum(is.na(Yt)), "\n")
cat("Final deltat length:", length(deltat), "\n")

# The bug is that when sum(deltat == 0) >= 10, the cleanup doesn't happen
# but downstream TMLE code expects clean data
if (sum(is.na(Yt)) > 0) {
  cat("BUG REPRODUCED: Missing values remain in Yt, which will cause TMLE to fail\n")
} else {
  cat("No missing values in Yt - this case works\n")
}

cat("\nSUMMARY:\n")
cat("- When there are < 10 missing values: cleanup runs, TMLE gets clean data\n")
cat("- When there are >= 10 missing values: cleanup doesn't run, TMLE gets dirty data and fails\n")
cat("- The fix is to remove the entire conditional block (lines 181-186 in vim-factors.R)\n")