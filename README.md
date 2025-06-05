
<!-- README.md is generated from README.Rmd. Please edit that file -->

# varimpact - variable importance through causal inference

[![Build
Status](https://travis-ci.org/ck37/varimpact.svg?branch=master)](https://travis-ci.org/ck37/varimpact)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ck37/varimpact?branch=master&svg=true)](https://ci.appveyor.com/project/ck37/varimpact)
[![codecov](https://codecov.io/gh/ck37/varimpact/branch/master/graph/badge.svg)](https://codecov.io/gh/ck37/varimpact)

## Summary

varimpact uses causal inference statistics to generate variable
importance estimates for a given dataset and outcome. It answers the
question: which of my Xs are most related to my Y? Each variable’s
influence on the outcome is estimated semiparametrically, without
assuming a linear relationship or other functional form, and the
covariate list is ranked by order of importance. This can be used for
exploratory data analysis, for dimensionality reduction, for
experimental design (e.g. to determine blocking and re-randomization),
to reduce variance in an estimation procedure, etc. See Hubbard,
Kennedy, & van der Laan (2018) for more details, or Hubbard & van der
Laan (2016) for an earlier description.

## Details

Each covariate is analyzed using targeted minimum loss-based estimation
([TMLE](https://CRAN.R-project.org/package=tmle)) as though it were a
treatment, with all other variables serving as adjustment variables via
[SuperLearner](https://github.com/ecpolley/SuperLearner). Then the
statistical significance of the estimated treatment effect for each
covariate determines the variable importance ranking. This formulation
allows the asymptotics of TMLE to provide valid standard errors and
p-values, unlike other variable importance algorithms.

The results provide raw p-values as well as p-values adjusted for false
discovery rate using the Benjamini-Hochberg (1995) procedure. Adjustment
variables are automatically clustered hierarchically using HOPACH (van
der Laan & Pollard 2003) in order to reduce dimensionality. The package
supports multi-core and multi-node parallelization, which are detected
and used automatically when a parallel backend is registered. Missing
values are automatically imputed using K-nearest neighbors (Troyanskaya
et al. 2001, Jerez et al. 2010) and missingness indicator variables are
incorporated into the analysis.

varimpact is under active development so please submit any bug reports
or feature requests to the [issue
queue](https://github.com/ck37/varimpact/issues), or email Alan and/or
Chris directly.

## Installation

### GitHub

``` r
# Install remotes if necessary:
# install.packages("remotes")
remotes::install_github("ck37/varimpact")
```

### CRAN

Forthcoming fall 2022

## Examples

### Example: basic functionality

``` r
library(varimpact)

####################################
# Create test dataset.
set.seed(1, "L'Ecuyer-CMRG")
N <- 200
num_normal <- 4
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
# Add some missing data to X so we can test imputation.
for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA

####################################
# Basic example
vim <- varimpact(Y = Y, data = X)
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 4 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 4 numerics.

# Review consistent and significant results.
vim
#> Significant and consistent results:
#>       Type  Estimate            CI95      P-value Adj. p-value  Est. RR
#> V3 Ordered 0.4986136 (0.255 - 0.742) 2.935984e-05 0.0001174394 2.908162
#>         CI95 RR   P-value RR Adj. p-value RR
#> V3 (1.92 - 4.4) 2.214955e-07    8.859819e-07

# Look at all results.
vim$results_all
#>       Type    Estimate             CI95      P-value Adj. p-value   Est. RR
#> V3 Ordered  0.49861358  (0.255 - 0.742) 2.935984e-05 0.0001174394 2.9081617
#> V4 Ordered  0.21853793 (-0.167 - 0.604) 1.334006e-01 0.2668012809 1.4110231
#> V2 Ordered  0.04733746 (-0.276 - 0.371) 3.872064e-01 0.5162752138 1.0698709
#> V1 Ordered -0.10939221 (-0.494 - 0.275) 7.116162e-01 0.7116161584 0.8266168
#>           CI95 RR   P-value RR Adj. p-value RR Consistent
#> V3   (1.92 - 4.4) 2.214955e-07    8.859819e-07       TRUE
#> V4 (0.833 - 2.39) 1.000755e-01    2.001509e-01       TRUE
#> V2 (0.664 - 1.72) 3.905406e-01    5.207208e-01       TRUE
#> V1 (0.441 - 1.55) 7.238555e-01    7.238555e-01       TRUE

# Plot the V2 impact.
plot_var("V2", vim)
```

![](images/README-example_1-1.png)<!-- -->

``` r

# Generate latex tables with results.
exportLatex(vim)

# Clean up LaTeX files
cleanup_latex_files()
```

### Example: customize outcome and propensity score estimation

``` r
Q_lib = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.rpartPrune")
g_lib = c("SL.mean", "SL.glmnet")
set.seed(1, "L'Ecuyer-CMRG")
(vim = varimpact(Y = Y, data = X, Q.library = Q_lib, g.library = g_lib))
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 4 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 4 numerics.
#> Significant and consistent results:
#>       Type Estimate            CI95      P-value Adj. p-value  Est. RR
#> V3 Ordered  0.56401 (0.326 - 0.802) 1.749015e-06 6.996059e-06 3.644982
#>          CI95 RR   P-value RR Adj. p-value RR
#> V3 (2.34 - 5.69) 6.234554e-09    2.493822e-08
```

### Example: parallel via multicore

``` r
library(future)
plan("multisession")
vim = varimpact(Y = Y, data = X)
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 4 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 4 numerics.
```

### Example: mlbench breast cancer

``` r
data(BreastCancer, package = "mlbench")
data = BreastCancer

# Create a numeric outcome variable.
data$Y = as.integer(data$Class == "malignant")

# Use multicore parallelization to speed up processing.
plan("multisession")
(vim = varimpact(Y = data$Y, data = subset(data, select = -c(Y, Class, Id))))
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 9 
#> - Numeric variables: 0 
#> 
#> Estimating variable importance for 9 factors.
#> Significant and consistent results:
#>                Type  Estimate            CI95      P-value Adj. p-value
#> Bare.nuclei  Factor 0.6174459   (0.5 - 0.735) 0.000000e+00 0.000000e+00
#> Mitoses      Factor 0.4092028 (0.333 - 0.486) 0.000000e+00 0.000000e+00
#> Cl.thickness Factor 0.5245860 (0.382 - 0.667) 3.027578e-13 9.082735e-13
#> Cell.size    Factor 0.5650275 (0.395 - 0.735) 3.313050e-11 5.963490e-11
#>               Est. RR       CI95 RR   P-value RR Adj. p-value RR
#> Bare.nuclei  3.682218 (2.21 - 6.14) 0.000000e+00    0.000000e+00
#> Mitoses      2.093929 (1.85 - 2.37) 3.023193e-11    1.360437e-10
#> Cl.thickness 2.952087 (2.13 - 4.08) 2.956850e-07    8.870549e-07
#> Cell.size    3.445132    (1.98 - 6) 4.465977e-06    8.038759e-06
plot_var("Mitoses", vim)
```

![](images/README-example_5-1.png)<!-- -->

## Authors

Alan E. Hubbard and Chris J. Kennedy, University of California, Berkeley

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. Journal of
the royal statistical society. Series B (Methodological), 289-300.

Gruber, S., & van der Laan, M. J. (2012). tmle: An R Package for
Targeted Maximum Likelihood Estimation. Journal of Statistical Software,
51(i13).

Hubbard, A. E., Kennedy, C. J., van der Laan, M. J. (2018).
Data-adaptive target parameters. In M. van der Laan & S. Rose (2018)
Targeted Learning in Data Science. Springer.

Hubbard, A. E., Kherad-Pajouh, S., & van der Laan, M. J. (2016).
Statistical Inference for Data Adaptive Target Parameters. The
international journal of biostatistics, 12(1), 3-19.

Hubbard, A., Munoz, I. D., Decker, A., Holcomb, J. B., Schreiber, M. A.,
Bulger, E. M., … & Rahbar, M. H. (2013). Time-Dependent Prediction and
Evaluation of Variable Importance Using SuperLearning in High
Dimensional Clinical Data. The journal of trauma and acute care surgery,
75(1 0 1), S53.

Hubbard, A. E., & van der Laan, M. J. (2016). Mining with inference:
data-adaptive target parameters (pp. 439-452). In P. Bühlmann et
al. (Ed.), Handbook of Big Data. CRC Press, Taylor & Francis Group, LLC:
Boca Raton, FL.

Jerez, J. M., Molina, I., García-Laencina, P. J., Alba, E., Ribelles,
N., Martín, M., & Franco, L. (2010). Missing data imputation using
statistical and machine learning methods in a real breast cancer
problem. Artificial intelligence in medicine, 50(2), 105-115.

Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular
and irregular histograms by penalized likelihood. Computational
Statistics & Data Analysis, 54(12), 3313-3323.

Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T.,
Tibshirani, R., Botstein, D., & Altman, R. B. (2001). Missing value
estimation methods for DNA microarrays. Bioinformatics, 17(6), 520-525.

van der Laan, M. J. (2006). Statistical inference for variable
importance. The International Journal of Biostatistics, 2(1).

van der Laan, M. J., & Pollard, K. S. (2003). A new algorithm for hybrid
hierarchical clustering with visualization and the bootstrap. Journal of
Statistical Planning and Inference, 117(2), 275-303.

van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super
learner. Statistical applications in genetics and molecular biology,
6(1).

van der Laan, M. J., & Rose, S. (2011). Targeted learning: causal
inference for observational and experimental data. Springer Science &
Business Media.
