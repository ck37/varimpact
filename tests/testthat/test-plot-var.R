library(varimpact)
library(testthat)

# plot_var() is exported and, until now, had no test at all. These fit one
# small varimpact model with a numeric and a factor variable and check the
# plot it produces for each, plus the two ways it refuses.

context("plot_var")

future::plan("sequential")

set.seed(3, "L'Ecuyer-CMRG")
N = 120
X = data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               f1 = factor(sample(c("a", "b", "c"), N, replace = TRUE)))
Y = rbinom(N, 1, plogis(0.6 * X$x1 - 0.3 * (X$f1 == "c")))

vim = varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
                Q.library = c("SL.mean", "SL.glm"),
                g.library = c("SL.mean", "SL.glm"),
                bins_numeric = 3L)

# The bar layer is the first geom; ggplot_build() is the stable way to read
# what was drawn without depending on ggplot2's internal object layout.
bar_data = function(p) ggplot2::ggplot_build(p)$data[[1]]

test_that("a numeric variable produces one bar per level plus the impact bar", {
  skip_if(is.null(vim$numeric_vims$results_by_level), "no numeric results")
  var = vim$numeric_vims$results_by_level$name[1]
  skip_if(!var %in% rownames(vim$results_all), "variable not in results_all")

  p = plot_var(var, vim)
  expect_s3_class(p, "ggplot")

  n_levels = sum(vim$numeric_vims$results_by_level$name == var)
  bars = bar_data(p)
  expect_equal(nrow(bars), n_levels + 1L)

  # The impact bar carries the overall estimate for the variable.
  estimate = vim$results_all[var, "Estimate"]
  expect_true(any(abs(bars$y - estimate) < 1e-8))

  # Low risk, high risk and impact are each drawn in their own colour.
  expect_gte(length(unique(bars$fill)), 3L)
})

test_that("a factor variable can be plotted too", {
  skip_if(is.null(vim$factor_vims$results_by_level), "no factor results")
  var = vim$factor_vims$results_by_level$name[1]
  skip_if(!var %in% rownames(vim$results_all), "variable not in results_all")

  p = plot_var(var, vim)
  expect_s3_class(p, "ggplot")
  n_levels = sum(vim$factor_vims$results_by_level$name == var)
  expect_equal(nrow(bar_data(p)), n_levels + 1L)
})

test_that("digits controls the rounding of the bar labels", {
  skip_if(is.null(vim$numeric_vims$results_by_level), "no numeric results")
  var = vim$numeric_vims$results_by_level$name[1]
  skip_if(!var %in% rownames(vim$results_all), "variable not in results_all")

  labels = ggplot2::ggplot_build(plot_var(var, vim, digits = 1L))$data[[2]]$label
  expect_true(all(labels == round(labels, 1L)))
})

test_that("an unknown variable is refused by name", {
  expect_error(plot_var("no_such_variable", vim), "no variable called")
})

test_that("a result object without per-level results is refused", {
  empty = list(numeric_vims = list(results_by_level = NULL),
               factor_vims = list(results_by_level = NULL))
  expect_error(plot_var("x1", empty), "No results_by_level")
})
