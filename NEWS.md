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
  (`tmle::tmle()` applies it once). How badly this behaved depended on the range
  of the outcome, because `plogis()` only saturates once its argument is far
  from zero:

  - An outcome far from zero, such as one ranging 31 to 68, saturated
    completely. Every prediction collapsed to the upper bound, which made the
    training estimate identical in every bin, so `which.max()` and `which.min()`
    selected the same bin, every fold was discarded as "min and max level are
    the same", and `varimpact()` returned no results at all.
  - An outcome near zero, such as one ranging -2 to 3, did not collapse. It
    came back distorted instead: plausible-looking estimates that were simply
    wrong, which is the harder case to notice.

  For a binary outcome the bounds are `c(0, 1)` and the extra transform was
  `plogis()` of a probability, which shifted theta into (0.5, 0.731) but left
  bin selection intact.

## Dependency changes

* Removed the `multtest`, `cvTools`, `caret` and `modeest` imports, and moved
  `glmnet` and `MASS` to Suggests (issue #65, first part). This drops one of
  the two Bioconductor dependencies and shrinks the recursive dependency tree
  from 117 packages to 69. What each import did is now done in a few lines of
  base R:

  - Holm and Benjamini-Hochberg adjusted p-values come from
    `stats::p.adjust()`, in a new internal `adjust_pvalues()`. The values are
    identical to `multtest::mt.rawp2adjp()`'s, including its treatment of a
    missing p-value as still counting toward the number of tests.
  - Median and knn imputation are done by the new internal `impute_median()`
    and `impute_knn()`, the latter calling `RANN::nn2()` directly, which
    `caret` did underneath. Both reproduce `caret::preProcess()`'s output to
    the bit. `impute_knn()` no longer errors when fewer than five complete
    rows exist; it uses as many as there are.
  - `modeest` was only referenced from a commented-out line.
  - `glmnet` and `MASS` were never called from this package; they are
    dependencies of SuperLearner wrappers a user may choose, and `SL.glmnet`
    is used in the tests, hence Suggests.

* Cross-validation folds are now genuinely random. `create_cv_folds()` used
  `cvTools::cvFolds(type = "random")$which`, but `$which` is the fold of the
  *permuted* observation, `rep(seq_len(V), length.out = n)`; the permutation
  itself sits in `$subsets`, which was never read. So every run assigned
  observations to folds by their row order within each outcome stratum,
  1, 2, ..., V, 1, 2, ... regardless of the seed. Folds are now drawn with
  `sample.int()`, still stratified on a binary outcome and still balanced to
  within one observation. Results of any given run therefore differ from
  those of earlier versions, and now depend on the seed as they were always
  meant to. A stratum with fewer observations than folds no longer errors.

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
