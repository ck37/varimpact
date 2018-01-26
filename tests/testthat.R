library(varimpact)

# Only run tests if testthat package is installed.
# This is in compliance with "Writing R Extensions" §1.1.3.1.
if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_check("varimpact", reporter = "check")
}
