# Solution for Issue #8: Transform continuous outcome estimates back to original scale

## Problem Summary
Issue #8 reported that continuous outcome estimates were being reported on the [0, 1] scale rather than the original scale of the outcome variable. This made the results difficult to interpret.

## Root Cause Analysis
The issue occurs because:
1. For continuous outcomes, TMLE transforms Y to [0,1] scale using: `Y_star = (Y - Qbounds[1]) / diff(Qbounds)`
2. All estimation is done on this transformed scale
3. Final estimates were not being transformed back to the original scale

## Solution Implemented

### 1. Modified `R/estimate_pooled_results.R`

**Added new parameters:**
```r
estimate_pooled_results = function(fold_results, verbose = FALSE, 
                                   Qbounds = NULL, map_to_ystar = FALSE)
```

**Added transformation logic:**
```r
# Transform thetas back to original scale if needed
if (map_to_ystar && !is.null(Qbounds) && length(Qbounds) == 2) {
  if (verbose) {
    cat("Transforming estimates from [0,1] scale back to original scale.\n")
    cat("Original Qbounds:", Qbounds, "\n")
    cat("Estimates before transformation:", range(thetas, na.rm = TRUE), "\n")
  }
  
  # Transform back: Y = Y_star * diff(Qbounds) + Qbounds[1]
  thetas = thetas * diff(Qbounds) + Qbounds[1]
  
  if (verbose) {
    cat("Estimates after transformation:", range(thetas, na.rm = TRUE), "\n")
  }
}
```

**Also transforms influence curve components:**
```r
# Also transform Q_star and Y_star back to original scale for influence curves
if (map_to_ystar && !is.null(Qbounds) && length(Qbounds) == 2) {
  combined_df$Q_star = combined_df$Q_star * diff(Qbounds) + Qbounds[1]
  combined_df$Y_star = combined_df$Y_star * diff(Qbounds) + Qbounds[1]
}
```

### 2. Modified `R/vim-numerics.R`

**Added logic to determine when transformation is needed:**
```r
# Determine if we need to transform back to original scale
# For continuous outcomes, the estimates are on [0,1] scale and need to be transformed back
map_to_ystar = FALSE
if (!is.null(Qbounds) && length(Qbounds) == 2) {
  # Check if this is a continuous outcome (not binary)
  if (family == "gaussian" || (family == "binomial" && length(unique(Y)) > 2)) {
    map_to_ystar = TRUE
  }
}
```

**Updated function calls:**
```r
pooled_min = estimate_pooled_results(lapply(fold_results, function(x) x$level_min),
                                     verbose = verbose,
                                     Qbounds = Qbounds,
                                     map_to_ystar = map_to_ystar)
pooled_max = estimate_pooled_results(lapply(fold_results, function(x) x$level_max),
                                     verbose = verbose,
                                     Qbounds = Qbounds,
                                     map_to_ystar = map_to_ystar)

pooled_bin = estimate_pooled_results(bin_list, verbose = verbose,
                                     Qbounds = Qbounds,
                                     map_to_ystar = map_to_ystar)
```

### 3. Modified `R/vim-factors.R`

**Applied the same changes as in vim-numerics.R:**
- Added logic to determine when transformation is needed
- Updated all calls to `estimate_pooled_results()` to pass the new parameters

## Testing

Created `test_continuous_outcome_fix.R` to verify the fix:
- Creates continuous outcome with known range (10-50)
- Runs varimpact with gaussian family
- Checks that estimates are on original scale, not [0,1] scale

## Files Changed

1. **R/estimate_pooled_results.R** - Core transformation logic
2. **R/vim-numerics.R** - Updated function calls for numeric variables
3. **R/vim-factors.R** - Updated function calls for factor variables
4. **test_continuous_outcome_fix.R** - Test script
5. **demonstrate_issue.md** - Documentation

## Expected Results

**Before fix:** Continuous outcome estimates between 0 and 1
**After fix:** Continuous outcome estimates on original scale (e.g., 10-50 range)

## Verification Steps

1. Create continuous outcome with known range
2. Run `varimpact(Y = Y_continuous, data = X, family = "gaussian")`
3. Check `vim$results_all$Estimate` values are on original scale
4. Check fold-specific estimates are also on original scale

This fix ensures that users get interpretable results on the original scale of their outcome variable, addressing the core issue reported in Issue #8.