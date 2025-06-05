# Fix for Issue #8: Transform continuous outcome estimates back to original scale

## Problem Description

The issue was that for continuous outcomes, varimpact was reporting estimates on the [0, 1] scale rather than the original scale of the outcome variable. This happened because:

1. For continuous outcomes, the TMLE algorithm transforms the outcome Y to the [0, 1] scale using the formula:
   ```
   Y_star = (Y - Qbounds[1]) / diff(Qbounds)
   ```
   where `Qbounds` are the bounds of the original outcome (extended by 10%).

2. All TMLE estimation is done on this [0, 1] scale.

3. The final estimates were not being transformed back to the original scale.

## Solution

The fix involves modifying the `estimate_pooled_results()` function to:

1. Accept additional parameters `Qbounds` and `map_to_ystar` to know when and how to transform back.

2. Transform the final theta estimates back to the original scale using:
   ```
   thetas_original = thetas * diff(Qbounds) + Qbounds[1]
   ```

3. Also transform the Q_star and Y_star values used in influence curve calculations back to the original scale.

## Files Modified

1. **R/estimate_pooled_results.R**: 
   - Added `Qbounds` and `map_to_ystar` parameters
   - Added transformation of thetas back to original scale
   - Added transformation of Q_star and Y_star for influence curve calculation

2. **R/vim-numerics.R**:
   - Updated calls to `estimate_pooled_results()` to pass bounds information
   - Added logic to determine when transformation is needed

3. **R/vim-factors.R**:
   - Updated calls to `estimate_pooled_results()` to pass bounds information
   - Added logic to determine when transformation is needed

## Test Case

The test case creates a continuous outcome with values roughly between 10 and 50, then checks that the final estimates are on this original scale rather than the [0, 1] scale.

Before the fix: Estimates would be between 0 and 1
After the fix: Estimates should be on the original scale (10-50 range)

## Verification

To verify the fix works:

1. Create a continuous outcome with a known range (e.g., 10-50)
2. Run varimpact with `family = "gaussian"`
3. Check that the estimates in `results_all$Estimate` are on the original scale
4. Check that fold-specific estimates are also on the original scale

The fix ensures that users get interpretable results on the original scale of their outcome variable, which is much more useful for practical applications.