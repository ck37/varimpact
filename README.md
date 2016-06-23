# varImpact
varImpact uses causal inference statistics to generate variable importance estimates for a given dataset and outcome. Each covariate is analyzed using targeted minimum loss-based estimation (TMLE) as though it were a treatment, with all other variables serving as adjustment variables via SuperLearner. Then the change in the outcome variable due to treatment is how the variable importance is estimated. This formulation allows the asymptotics of TMLE to provide valid standard errors and p-values, unlike other variable importance estimates. See Hubbard & van der Laan (2016) for more details.

The results provide raw p-values as well as p-values adjusted for false discovery rate using the Benjamini-Hochman (1995) procedure. Adjustment variables are automatically clustered hierarchically using HOPACH in order to reduce dimensionality.  The package supports multi-core and multi-node parallelization, which are detected and used automatically when a parallel backend is registered. Missing values are automatically imputed using K-nearest neighbors and missingness indicator variables are incorporated into the analysis.

varImpact is under active development so please submit any bug reports or feature requests to the issue queue, or email Alan & Chris directly.

## Install

### Bioconductor packages

A few package requirements are installed via [bioconductor](https://www.bioconductor.org) rather than CRAN.

First, [install bioconductor](https://www.bioconductor.org/install/) if you don't have it already:
```{r}
source("http://bioconductor.org/biocLite.R")
biocLite()
```

Then install the packages:
```{r}
biocLite(c("hopach", "multtest"))
```

### varImpact install

```{r}
# Install devtools if necessary:
if (!require("devtools")) install.packages("devtools")
devtools::install_github("ck37/varImpact")
```

### Examples

```{r}
# Create test dataset.
set.seed(1)
N = 200
num_normal = 7
X = as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
Y = rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
# Add some missing data to X.
miss_num = 10
for (i in 1:miss_num) X[sample(nrow(X), 1), sample(ncol(X), 1)] = NA

# Basic example
vim <- varImpact(Y = Y, data = X)
vim
vim$results_all
exportLatex(vim)

# Impute by median rather than knn.
vim <- varImpact(Y = Y, data = X, impute = "median")

# doMC parallel (multicore) example.
library(doMC)
registerDoMC()
vim <- varImpact(Y = Y, data = X)

# doSNOW parallel example.
library(doSNOW)
library(RhpcBLASctl)
# Detect the number of physical cores on this computer using RhpcBLASctl.
cluster <- makeCluster(get_num_cores())
registerDoSNOW(cluster)
vim <- varImpact(Y = Y, data = X)
stopCluster(cluster)
```

## Authors

Alan Hubbard and Chris Kennedy, University of California, Berkeley

## References

Gruber, S., & Laan, M. V. D. (2012). tmle: An R Package for Targeted Maximum Likelihood Estimation. Journal of Statistical Software, 51(i13).

Hubbard, A., & van der Laan, M. (2016). Mining with inference: data-adaptive target parameter (pp. 439-452). In P. Bühlmann et al. (Ed.), Handbook of Big Data. CRC Press, Taylor & Francis Group, LLC: Boca Raton, FL.

Van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. Statistical applications in genetics and molecular biology, 6(1).

Van der Laan, M. J., & Rose, S. (2011). Targeted learning: causal inference for observational and experimental data. Springer Science & Business Media.
