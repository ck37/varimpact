# varimpact 1.3.0-9007 (development version)

## Bug fixes

* Fixed issue #8: variable importance estimates for a continuous outcome are now
  reported on the scale of the outcome rather than the internal [0, 1] scale.
  `varimpact()` maps a continuous `Y` into [0, 1] with
  `Y_star = (Y - Qbounds[1]) / diff(Qbounds)` before running the CV-TMLE;
  `estimate_pooled_results()` now applies the inverse of that map to the
  fluctuated `Q_star`, so the treatment-specific means, the risk difference, the
  risk ratio and the influence curves all come back on the original scale. For a
  binary outcome `Qbounds` is `c(0, 1)` and the transformation is the identity,
  so binary results are unchanged.

* Fixed a duplicated line in `estimate_tmle2()` that applied `plogis()` to
  `Qstar` twice when mapping the targeted predictions back to the outcome scale
  (`tmle::tmle()` applies it once). For a continuous outcome the second
  `plogis()` saturated every prediction, which made the training estimate
  identical in every bin; `which.max()` and `which.min()` then selected the same
  bin, every fold was discarded as "min and max level are the same", and
  `varimpact()` returned no results at all for `family = "gaussian"`.

## Internal changes

* `estimate_pooled_results()` gained a `Qbounds` argument, defaulting to
  `c(0, 1)`. `vim_numerics()` and `vim_factors()` pass down the bounds that
  `varimpact()` computed.

* Added `tests/testthat/test-continuous-outcome-scale.R`, covering the exact
  rescaling of thetas and influence curves, an end-to-end gaussian run, and a
  binary run as a regression check.

## Known limitations

* The risk ratio (`EY1 / EY0`) is only interpretable for an outcome that is
  strictly positive. Now that continuous estimates are on the original scale, an
  outcome whose range straddles zero can produce a negative or undefined risk
  ratio; the risk difference columns are unaffected.

# varimpact 1.2.4

Previous releases...
