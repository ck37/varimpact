
<!-- README.md is generated from README.Rmd. Please edit that file -->
varImpact - variable importance through causal inference
========================================================

[![Build Status](https://travis-ci.org/ck37/varImpact.svg?branch=master)](https://travis-ci.org/ck37/varImpact) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ck37/varImpact?branch=master&svg=true)](https://ci.appveyor.com/project/ck37/varImpact) [![codecov](https://codecov.io/gh/ck37/varImpact/branch/master/graph/badge.svg)](https://codecov.io/gh/ck37/varImpact)

Summary
-------

varImpact uses causal inference statistics to generate variable importance estimates for a given dataset and outcome. It answers the question: which of my Xs are most related to my Y? Each variable's influence on the outcome is estimated semiparametrically, without assuming a linear relationship or other functional form, and the covariate list is ranked by order of importance. This can be used for exploratory data analysis, for dimensionality reduction, for experimental design (e.g. to determine blocking and re-randomization), to reduce variance in an estimation procedure, etc. See Hubbard & van der Laan (2016) and Hubbard, Kennedy, & van der Laan (2017) for more details.

Details
-------

Each covariate is analyzed using targeted minimum loss-based estimation ([TMLE](https://cran.r-project.org/web/packages/tmle/index.html)) as though it were a treatment, with all other variables serving as adjustment variables via [SuperLearner](https://github.com/ecpolley/SuperLearner). Then the statistical significance of the estimated treatment effect for each covariate determines the variable importance ranking. This formulation allows the asymptotics of TMLE to provide valid standard errors and p-values, unlike other variable importance algorithms.

The results provide raw p-values as well as p-values adjusted for false discovery rate using the Benjamini-Hochberg (1995) procedure. Adjustment variables are automatically clustered hierarchically using HOPACH (van der Laan & Pollard 2003) in order to reduce dimensionality. The package supports multi-core and multi-node parallelization, which are detected and used automatically when a parallel backend is registered. Missing values are automatically imputed using K-nearest neighbors (Troyanskaya et al. 2001, Jerez et al. 2010) and missingness indicator variables are incorporated into the analysis.

varImpact is under active development so please submit any bug reports or feature requests to the [issue queue](https://github.com/ck37/varImpact/issues), or email Alan & Chris directly.

Installation
------------

### Github

``` r
# Install devtools if necessary:
# install.packages("devtools")
devtools::install_github("ck37/varImpact")
```

### CRAN

Forthcoming Fall 2017

Examples
--------

``` r
library(varImpact)
#> Loading required package: SuperLearner
#> Loading required package: nnls
#> Super Learner
#> Version: 2.0-23-9000
#> Package created on 2017-07-20

####################################
# Create test dataset.
set.seed(1)
N <- 200
num_normal <- 5
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
# Add some missing data to X so we can test imputation.
for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA

####################################
# Basic example
vim <- varImpact(Y = Y, data = X)
#> Finished pre-processing variables.
#> 
#> Processing results:
#> - Factor variables: 0 
#> - Numeric variables: 5 
#> 
#> No factor variables - skip VIM estimation.
#> 
#> Estimating variable importance for 5 numerics.
#> Loading required package: glmnet
#> Loading required package: Matrix
#> Loading required package: foreach
#> Loaded glmnet 2.0-13
vim
#> No significant and consistent results.
#> All results:
#>       Type   Estimate              CI95    P-value Adj. p-value Consistent
#> V2 Ordered 0.16418106 (-0.0482 - 0.377) 0.06483449    0.1783680       TRUE
#> V4 Ordered 0.15611465 (-0.0712 - 0.383) 0.08918400    0.1783680       TRUE
#> V3 Ordered 0.10041232  (-0.133 - 0.334) 0.19933203    0.2657760       TRUE
#> V5 Ordered 0.04340259  (-0.196 - 0.283) 0.36142148    0.3614215       TRUE
vim$results_all
#>       Type   Estimate              CI95    P-value Adj. p-value Consistent
#> V2 Ordered 0.16418106 (-0.0482 - 0.377) 0.06483449    0.1783680       TRUE
#> V4 Ordered 0.15611465 (-0.0712 - 0.383) 0.08918400    0.1783680       TRUE
#> V3 Ordered 0.10041232  (-0.133 - 0.334) 0.19933203    0.2657760       TRUE
#> V5 Ordered 0.04340259  (-0.196 - 0.283) 0.36142148    0.3614215       TRUE
exportLatex(vim)
#> NULL
# Clean up - will get a warning if there were no consistent results.
suppressWarnings({
  file.remove(c("varimpByFold.tex", "varImpAll.tex", "varimpConsistent.tex"))
})
#> [1]  TRUE  TRUE FALSE

# Impute by median rather than knn.
vim <- varImpact(Y = Y, data = X, impute = "median")
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
vim <- varImpact(Y = Y, data = X, Q.library = Q_lib, g.library = g_lib)
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
#> No significant and consistent results.
#> All results:
#>       Type  Estimate              CI95   P-value Adj. p-value Consistent
#> V2 Ordered 0.1626026 (-0.0505 - 0.376) 0.0674013    0.2478408       TRUE
#> V4 Ordered 0.1419863 (-0.0988 - 0.383) 0.1239204    0.2478408       TRUE
#> V3 Ordered 0.1011522  (-0.132 - 0.334) 0.1971035    0.2628047       TRUE
#> V5 Ordered 0.0417297  (-0.197 - 0.281) 0.3661919    0.3661919       TRUE

####################################
# Parallel (multicore) example.
library(future)
plan("multiprocess")
vim <- varImpact(Y = Y, data = X)
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
vim <- varImpact(Y = Y, data = X)
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
vim <- varImpact(Y = data$Y, data = subset(data, select = -c(Y, Class, Id)))
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
#> Mitoses       Factor 0.2656178 1.776357e-15 1.598721e-14 (0.199 - 0.332)
#> Marg.adhesion Factor 0.3944912 7.468470e-13 3.360812e-12 (0.285 - 0.504)
#> Cl.thickness  Factor 0.3775814 3.032512e-08 6.823152e-08 (0.241 - 0.514)
#> Cell.size     Factor 0.2702853 1.972746e-05 2.536387e-05 (0.141 - 0.399)
```

Authors
-------

Alan E. Hubbard and Chris J. Kennedy, University of California, Berkeley

References
----------

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.

Gruber, S., & van der Laan, M. J. (2012). tmle: An R Package for Targeted Maximum Likelihood Estimation. Journal of Statistical Software, 51(i13).

Hubbard, A. E., Kennedy, C. J., van der Laan, M. J. (2017). Data-adaptive Variable Importance Using Cross-validated Targeted Maximum Likelihood Estimation. In M. van der Laan & S. Rose (2017) Targeted Learning in Data Science. Springer.

Hubbard, A. E., Kherad-Pajouh, S., & van der Laan, M. J. (2016). Statistical Inference for Data Adaptive Target Parameters. The international journal of biostatistics, 12(1), 3-19.

Hubbard, A., Munoz, I. D., Decker, A., Holcomb, J. B., Schreiber, M. A., Bulger, E. M., ... & Rahbar, M. H. (2013). Time-Dependent Prediction and Evaluation of Variable Importance Using SuperLearning in High Dimensional Clinical Data. The journal of trauma and acute care surgery, 75(1 0 1), S53.

Hubbard, A. E., & van der Laan, M. J. (2016). Mining with inference: data-adaptive target parameters (pp. 439-452). In P. Bühlmann et al. (Ed.), Handbook of Big Data. CRC Press, Taylor & Francis Group, LLC: Boca Raton, FL.

Jerez, J. M., Molina, I., García-Laencina, P. J., Alba, E., Ribelles, N., Martín, M., & Franco, L. (2010). Missing data imputation using statistical and machine learning methods in a real breast cancer problem. Artificial intelligence in medicine, 50(2), 105-115.

Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. Computational Statistics & Data Analysis, 54(12), 3313-3323.

Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., Tibshirani, R., Botstein, D., & Altman, R. B. (2001). Missing value estimation methods for DNA microarrays. Bioinformatics, 17(6), 520-525.

van der Laan, M. J. (2006). Statistical inference for variable importance. The International Journal of Biostatistics, 2(1).

van der Laan, M. J., & Pollard, K. S. (2003). A new algorithm for hybrid hierarchical clustering with visualization and the bootstrap. Journal of Statistical Planning and Inference, 117(2), 275-303.

van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. Statistical applications in genetics and molecular biology, 6(1).

van der Laan, M. J., & Rose, S. (2011). Targeted learning: causal inference for observational and experimental data. Springer Science & Business Media.
