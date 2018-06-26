
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
to reduce variance in an estimation procedure, etc. See Hubbard & van
der Laan (2016) and Hubbard, Kennedy, & van der Laan (2017) for more
details.

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
queue](https://github.com/ck37/varimpact/issues), or email Alan & Chris
directly.

## Installation

### Github

``` r
# Install devtools if necessary:
# install.packages("devtools")
devtools::install_github("ck37/varimpact")
```

### CRAN

Forthcoming Summer 2018

## Examples

``` r
library(varimpact)
#> Loading required package: SuperLearner
#> Loading required package: nnls
#> Super Learner
#> Version: 2.0-24-9000
#> Package created on 2018-03-09

####################################
# Create test dataset.
set.seed(1, "L'Ecuyer-CMRG")
N <- 200
num_normal <- 5
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
#> - Numeric variables: 5 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 5 numerics.
vim
#> Significant and consistent results:
#>       Type  Estimate     P-value Adj. P-value            CI 95
#> V2 Ordered 0.2501553 0.004721456   0.02360728 (0.0613 - 0.439)
vim$results_all
#>       Type    Estimate              CI95     P-value Adj. p-value
#> V2 Ordered  0.25015526  (0.0613 - 0.439) 0.004721456   0.02360728
#> V5 Ordered  0.14126628 (-0.0724 - 0.355) 0.097520684   0.24380171
#> V3 Ordered  0.15309902  (-0.136 - 0.443) 0.150010821   0.25001804
#> V1 Ordered  0.03420408  (-0.224 - 0.293) 0.397651895   0.49706487
#> V4 Ordered -0.16142012 (-0.383 - 0.0602) 0.923254928   0.92325493
#>    Consistent
#> V2       TRUE
#> V5       TRUE
#> V3      FALSE
#> V1       TRUE
#> V4      FALSE
exportLatex(vim)
# Clean up - will get a warning if there were no consistent results.
suppressWarnings({
  file.remove(c("varimpByFold.tex", "varImpAll.tex", "varimpConsistent.tex"))
})
#> [1] TRUE TRUE TRUE

# Impute by median rather than knn.
vim <- varimpact(Y = Y, data = X, impute = "median")
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 5 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 5 numerics.

# Customize Q and g libraries for TMLE estimation.
Q_lib <- c("SL.mean", "SL.glmnet", "SL.ranger", "SL.rpartPrune", "SL.bayesglm")
g_lib <- c("SL.mean", "SL.glmnet")
vim <- varimpact(Y = Y, data = X, Q.library = Q_lib, g.library = g_lib)
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 5 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 5 numerics.
vim
#> Significant and consistent results:
#>       Type  Estimate     P-value Adj. P-value            CI 95
#> V2 Ordered 0.2563156 0.004167914   0.02083957 (0.0659 - 0.447)

####################################
# Parallel (multicore) example.
library(future)
plan("multiprocess")
vim <- varimpact(Y = Y, data = X)
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 5 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 5 numerics.

####################################
# SNOW parallel example.
library(RhpcBLASctl)
# Detect the number of physical cores on this computer using RhpcBLASctl.
cl = parallel::makeCluster(get_num_cores())
plan("cluster", workers = cl)
vim <- varimpact(Y = Y, data = X)
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 5 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 5 numerics.
parallel::stopCluster(cl)

####################################
# mlbench BreastCancer example.
data(BreastCancer, package = "mlbench")
data <- BreastCancer

# Create a numeric outcome variable.
data$Y <- as.numeric(data$Class == "malignant")

# Use multicore parallelization to speed up processing.
plan("multiprocess")
vim <- varimpact(Y = data$Y, data = subset(data, select = -c(Y, Class, Id)))
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 9 
#> - Numeric variables: 0 
#> 
#> Estimating variable importance for 9 factors.
#> No numeric variables for variable importance estimation.
vim
#> Significant and consistent results:
#>                 Type  Estimate      P-value Adj. P-value           CI 95
#> Mitoses       Factor 0.3177571 0.000000e+00 0.000000e+00  (0.255 - 0.38)
#> Marg.adhesion Factor 0.3807888 5.184742e-14 2.333134e-13  (0.28 - 0.481)
#> Cell.size     Factor 0.5512742 2.247580e-10 6.742740e-10 (0.378 - 0.725)
#> Cl.thickness  Factor 0.3793746 8.141966e-09 1.465554e-08 (0.248 - 0.511)
```

## Authors

Alan E. Hubbard and Chris J. Kennedy, University of California, Berkeley

## References

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. Journal of
the royal statistical society. Series B (Methodological), 289-300.

Gruber, S., & van der Laan, M. J. (2012). tmle: An R Package for
Targeted Maximum Likelihood Estimation. Journal of Statistical Software,
51(i13).

Hubbard, A. E., Kennedy, C. J., van der Laan, M. J. (2017).
Data-adaptive Variable Importance Using Cross-validated Targeted Maximum
Likelihood Estimation. In M. van der Laan & S. Rose (2017) Targeted
Learning in Data Science. Springer.

Hubbard, A. E., Kherad-Pajouh, S., & van der Laan, M. J. (2016).
Statistical Inference for Data Adaptive Target Parameters. The
international journal of biostatistics, 12(1), 3-19.

Hubbard, A., Munoz, I. D., Decker, A., Holcomb, J. B., Schreiber, M. A.,
Bulger, E. M., … & Rahbar, M. H. (2013). Time-Dependent Prediction and
Evaluation of Variable Importance Using SuperLearning in High
Dimensional Clinical Data. The journal of trauma and acute care surgery,
75(1 0 1), S53.

Hubbard, A. E., & van der Laan, M. J. (2016). Mining with inference:
data-adaptive target parameters (pp. 439-452). In P. Bühlmann et al.
(Ed.), Handbook of Big Data. CRC Press, Taylor & Francis Group, LLC:
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
