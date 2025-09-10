# Force the bug condition by creating a scenario with >= 10 missing values in training fold

set.seed(1, "L'Ecuyer-CMRG")
N <- 200
num_normal <- 4
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))

# Add missing data to Y - create many missing values to ensure >= 10 in training fold
missing_indices <- c(4,6,7,8,11,15,20,21,28,32,72,80,85,90,95,100,105,110,115,120)
Y[missing_indices] <- NA

cat("Total missing Y values:", sum(is.na(Y)), "\n")

# Use single fold to ensure all missing values are in training data
folds <- rep(1, N)  # All data in fold 1
folds[1:10] <- 2    # Only first 10 observations in fold 2
fold_k <- 2         # Use fold 2 as validation, so most data (including missing) is in training

# Simulate training data (all data not in fold 2, i.e., most of the data)
Yt <- Y[folds != fold_k]
At <- as.factor(sample(c("A", "B", "C"), sum(folds != fold_k), replace = TRUE))
Wtsht <- matrix(rnorm(sum(folds != fold_k) * 3), ncol = 3)

cat("Training fold size:", length(Yt), "\n")

# This is the key line from vim-factors.R line 175
deltat <- as.numeric(!is.na(Yt) & !is.na(At))

cat("Number of missing observations in training fold (deltat == 0):", sum(deltat == 0), "\n")

# Test both conditions
cat("\n=== TESTING THE BUG CONDITION ===\n")
cat("sum(deltat == 0) =", sum(deltat == 0), "\n")
cat("sum(deltat == 0) < 10 =", sum(deltat == 0) < 10, "\n")

if (sum(deltat == 0) < 10) {
  cat("CASE 1: CLEANUP CODE RUNS (< 10 missing)\n")
  cat("This is the working case - missing observations are removed\n")
} else {
  cat("CASE 2: CLEANUP CODE DOES NOT RUN (>= 10 missing)\n")
  cat("THIS IS THE BUG! Missing observations remain in data\n")
  cat("TMLE will receive data with missing values and fail\n")
}

# Show the problematic code behavior
cat("\n=== SIMULATING THE PROBLEMATIC CODE ===\n")
Yt_before <- Yt
At_before <- At
Wtsht_before <- Wtsht
deltat_before <- deltat

# This is the exact problematic code from vim-factors.R lines 181-186
if (sum(deltat == 0) < 10) {
  Yt <- Yt[deltat == 1]
  At <- At[deltat == 1]
  Wtsht <- Wtsht[deltat == 1, , drop = FALSE]
  deltat <- deltat[deltat == 1]
  cat("Cleanup performed: removed", sum(deltat_before == 0), "missing observations\n")
} else {
  cat("No cleanup performed: ", sum(deltat == 0), "missing observations remain\n")
}

cat("Before cleanup - Yt length:", length(Yt_before), "missing:", sum(is.na(Yt_before)), "\n")
cat("After cleanup  - Yt length:", length(Yt), "missing:", sum(is.na(Yt)), "\n")

if (sum(is.na(Yt)) > 0) {
  cat("\n*** BUG REPRODUCED ***\n")
  cat("Missing values remain in Yt, which will cause TMLE estimation to fail\n")
  cat("This happens when there are >= 10 missing observations\n")
} else {
  cat("\nNo bug in this case - all missing values were cleaned up\n")
}

cat("\n=== SOLUTION ===\n")
cat("Remove the entire conditional block (lines 181-186 in vim-factors.R):\n")
cat("if (sum(deltat == 0) < 10) {\n")
cat("  Yt = Yt[deltat == 1]\n")
cat("  At = At[deltat == 1]\n")
cat("  Wtsht = Wtsht[deltat == 1, , drop = FALSE]\n")
cat("  deltat = deltat[deltat == 1]\n")
cat("}\n")
cat("This will allow proper delta missingness estimation as mentioned in the TODO comment.\n")