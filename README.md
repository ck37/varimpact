# varImpact
Variable importance using causal inference

Author: Alan Hubbard

## Install

### Bioconductor packages

A few package requirements are installed via [bioconductor](https://www.bioconductor.org) rather than CRAN.

First, [install bioconductor](https://www.bioconductor.org/install/) if you don't have it already:
```{r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
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
