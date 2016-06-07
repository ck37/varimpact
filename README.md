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
