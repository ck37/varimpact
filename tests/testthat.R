library(testthat)
library(varImpact)
# We need to explictly load SuperLearner due to an issue with
# the "All" screener as of 2016-12-06.
library(SuperLearner)

test_check("varImpact")
