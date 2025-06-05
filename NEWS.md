# varimpact 1.2.5 (development version)

## Bug fixes

* Fixed Issue #8: Continuous outcome estimates are now correctly transformed back to the original scale instead of being reported on the [0, 1] scale. This makes results much more interpretable for users working with continuous outcomes. The fix modifies `estimate_pooled_results()` to accept bounds information and transform final estimates using the inverse of the TMLE transformation: `thetas = thetas * diff(Qbounds) + Qbounds[1]`.

## Internal changes

* Added `Qbounds` and `map_to_ystar` parameters to `estimate_pooled_results()` function
* Updated `vim_numerics()` and `vim_factors()` to pass bounds information for continuous outcome transformation
* Enhanced influence curve calculations to use original scale values

# varimpact 1.2.4

Previous releases...