# varimpact 1.3.0-9007 (development version)

## Bug fixes

* The relative-risk p-values (`P-value RR`, `Adj. p-value RR` in
  `results_all`; `rr_rawp`, `rr_Holm`, `rr_BH` in `results_raw`) were attached
  to the wrong variables whenever the relative-risk ranking differed from the
  risk-difference ranking. `compile_results()` ordered its rows by the
  risk-difference p-value but pasted in the relative-risk p-values as
  `multtest::mt.rawp2adjp()` returned them, sorted by their own order, so the
  `rr_rawp` column was always ascending regardless of which variable each row
  described. The relative-risk estimates and confidence intervals were placed
  correctly; only the three p-value columns were shifted. They now belong to
  their own variable. Found while replacing `multtest` (below); the README
  tables are re-rendered and show the corrected columns.

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
    missing p-value as still counting toward the number of tests. (The
    relative-risk p-value columns do change, because of the bug fix above.)
  - Median and knn imputation are done by the new internal `impute_median()`
    and `impute_knn()`, the latter calling `RANN::nn2()` directly, which
    `caret` did underneath. Both reproduce `caret::preProcess()`'s output to
    the bit. `impute_knn()` no longer errors when fewer than five complete
    rows exist; it uses as many as there are.
  - `modeest` was only referenced from a commented-out line.
  - `glmnet` and `MASS` were never called from this package; they are
    dependencies of SuperLearner wrappers a user may choose, and `SL.glmnet`
    is used in the tests, hence Suggests.

* Cross-validation fold assignment is unchanged, and that is worth spelling
  out. `create_cv_folds()` read `cvTools::cvFolds(type = "random")$which`, but
  `$which` is the fold of the *permuted* observation,
  `rep(seq_len(V), length.out = n)`; the permutation itself sits in
  `$subsets`, which was never read. So every run has always assigned
  observations to folds by their row order within each outcome stratum,
  1, 2, ..., V, 1, 2, ..., whatever the seed. The replacement reproduces that
  exactly. It also still makes (and discards) the same random draw per stratum
  that `cvFolds()` made, so that SuperLearner's internal folds and the
  per-variable seeds downstream see the same random numbers, and a fit from a
  given seed is identical to the previous version's. Making the folds
  genuinely random is left for a separate change: trying it on the mlbench
  BreastCancer example exposed a case where `hopach::hopach(mss = "mean")` in
  `reduce_dimensions()` never returns, so it needs a guard around HOPACH
  first. One difference: a stratum with fewer observations than folds no
  longer errors.

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
