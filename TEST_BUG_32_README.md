# Test for GitHub Issue #32: Missing Y Values Bug

## Bug Description

This test reproduces the bug reported in [GitHub issue #32](https://github.com/ck37/varimpact/issues/32).

### The Problem

The bug is located in `R/vim-factors.R` at lines 181-186:

```r
# TODO (CK): don't do this, in order to use the delta missingness estimation.
# To avoid crashing TMLE function just drop obs missing A or Y if the
# total number of missing is < 10
if (sum(deltat == 0) < 10) {
  Yt = Yt[deltat == 1]
  At = At[deltat == 1]
  Wtsht = Wtsht[deltat == 1, , drop = FALSE]
  deltat = deltat[deltat == 1]
}
```

### Root Cause

The problematic code creates inconsistent behavior:

1. **When there are < 10 missing values**: The cleanup code runs, removing missing observations before TMLE estimation
2. **When there are ≥ 10 missing values**: The cleanup code does NOT run, leaving missing values in the data
3. **Result**: TMLE estimation fails when there are ≥ 10 missing values because it receives uncleaned data

### The Fix

According to the TODO comment, this entire code block should be removed to allow proper delta missingness estimation.

## Test File

The test is located at: `tests/testthat/test-missing-y-bug.R`

### Test Cases

1. **Main Bug Test**: Tests with 11 missing Y values (should fail with current code)
2. **Working Case**: Tests with 3 missing Y values (works but uses problematic code path)  
3. **Edge Case**: Tests with exactly 10 missing Y values (should fail due to ≥ 10 condition)

## How to Run the Test

### Prerequisites

1. Install R and required dependencies:
```bash
# Install R
sudo apt-get install r-base

# Install required R packages
R -e "install.packages(c('testthat', 'SuperLearner', 'tmle', 'future', 'future.apply'), repos='https://cran.r-project.org')"
```

2. Install the varimpact package dependencies (this may take some time):
```r
# In R console
install.packages(c(
  'arules', 'caret', 'cvTools', 'dplyr', 'future', 'future.apply',
  'ggplot2', 'glmnet', 'histogram', 'hopach', 'magrittr', 'MASS',
  'modeest', 'multtest', 'RANN', 'tmle', 'xtable'
), repos='https://cran.r-project.org')
```

### Running the Test

```bash
cd /path/to/varimpact
R -e "library(testthat); test_file('tests/testthat/test-missing-y-bug.R')"
```

### Expected Results

With the current buggy code:
- ✅ Test with <10 missing Y values should PASS
- ❌ Test with exactly 10 missing Y values should FAIL  
- ❌ Test with >10 missing Y values should FAIL

After fixing the bug (removing lines 181-186 from vim-factors.R):
- ✅ All tests should PASS

## Reproduction Script

You can also run the bug reproduction directly:

```r
# Reproduce the exact simulation from the GitHub issue
set.seed(1, "L'Ecuyer-CMRG")
N <- 200
num_normal <- 4
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))

# Add some missing data to X
for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA

# Add missing data to Y - this triggers the bug (11 missing values)
Y[c(4,6,7,8,11,15,20,21,28,32,72)] <- NA

# This should fail with the current code
library(varimpact)
vim <- varimpact(Y = Y, data = X)
```

## Technical Details

### Why the Bug Occurs

1. `deltat = as.numeric(!is.na(Yt) & !is.na(At))` creates a vector where 1 = non-missing, 0 = missing
2. `sum(deltat == 0)` counts the number of missing observations
3. When this count is ≥ 10, the cleanup code is skipped
4. TMLE estimation receives data with missing values and fails

### The Inconsistency

The condition `sum(deltat == 0) < 10` creates an arbitrary threshold that leads to:
- Inconsistent data preprocessing
- Unpredictable failures based on the number of missing values
- Violation of the principle that similar inputs should produce similar behavior

### Proper Solution

Remove the entire conditional block (lines 181-186) to ensure consistent handling of missing data through the delta missingness estimation approach mentioned in the TODO comment.