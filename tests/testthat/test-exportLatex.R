# Test exportLatex functionality
library(testthat)
library(varimpact)

context("exportLatex() function")

# Skip exportLatex tests during R CMD check to avoid creating LaTeX files
skip_if(identical(Sys.getenv("_R_CHECK_PACKAGE_NAME_"), "varimpact"), 
        "Skipping exportLatex tests during R CMD check")

test_that("exportLatex creates and cleans up LaTeX files", {
  # Create test dataset
  set.seed(1, "L'Ecuyer-CMRG")
  N = 100
  X = data.frame(
    x1 = rnorm(N),
    x2 = rnorm(N),
    x3 = rnorm(N)
  )
  Y = rbinom(N, 1, plogis(0.2 * X$x1 + 0.1 * X$x2 - 0.2 * X$x3))
  
  # Run varimpact with minimal settings for speed
  future::plan("sequential")
  vim = varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
                  Q.library = c("SL.mean", "SL.glm"),
                  g.library = c("SL.mean", "SL.glm"),
                  bins_numeric = 3L)
  
  # Skip test if no results were generated (due to sample size constraints)
  skip_if(is.null(vim$results_all), "No varimpact results generated")
  
  # Define cleanup function to ensure files are always removed
  tex_files = c("varimpByFold.tex", "varimpAll.tex", "varimpConsistent.tex")
  cleanup_files = function() {
    cleanup_latex_files(verbose = FALSE)
  }
  
  # Ensure cleanup happens even if test fails
  on.exit(cleanup_files())
  
  # Test 1: exportLatex should create files
  exportLatex(vim)
  
  existing_files = tex_files[file.exists(tex_files)]
  
  expect_true(length(existing_files) > 0, 
              info = "exportLatex should create at least some LaTeX files")
  expect_true("varimpByFold.tex" %in% existing_files,
              info = "varimpByFold.tex should be created")
  expect_true("varimpAll.tex" %in% existing_files,
              info = "varimpAll.tex should be created")
  
  # Test 2: Manual cleanup should work
  cleanup_files()
  
  remaining_files = tex_files[file.exists(tex_files)]
  expect_equal(length(remaining_files), 0,
               info = "Manual cleanup should remove all LaTeX files")
  
  # Test 3: Manual cleanup after exportLatex should work
  exportLatex(vim)
  cleanup_files()
  
  remaining_files_after_cleanup = tex_files[file.exists(tex_files)]
  expect_equal(length(remaining_files_after_cleanup), 0,
               info = "Manual cleanup after exportLatex should remove LaTeX files")
})

test_that("exportLatex handles NULL results gracefully", {
  # Create a mock varimpact object with NULL results
  mock_vim = list(
    results_by_fold = NULL,
    results_all = NULL,
    results_consistent = data.frame()
  )
  
  # Should return NULL and give a warning
  expect_warning(
    result <- exportLatex(mock_vim),
    "Cannot export LaTeX: varimpact results are NULL or incomplete"
  )
  expect_null(result)
  
  # Should not create any files
  tex_files = c("varimpByFold.tex", "varimpAll.tex", "varimpConsistent.tex")
  existing_files = tex_files[file.exists(tex_files)]
  expect_equal(length(existing_files), 0,
               info = "exportLatex with NULL results should not create files")
})

test_that("exportLatex with custom outname and directory", {
  # Create test dataset
  set.seed(2, "L'Ecuyer-CMRG")
  N = 100
  X = data.frame(x1 = rnorm(N), x2 = rnorm(N))
  Y = rbinom(N, 1, plogis(0.3 * X$x1))
  
  # Run varimpact
  future::plan("sequential")
  vim = varimpact(Y = Y, data = X, V = 2L, verbose = FALSE,
                  Q.library = "SL.mean", g.library = "SL.mean",
                  bins_numeric = 3L)
  
  # Skip test if no results were generated
  skip_if(is.null(vim$results_all), "No varimpact results generated")
  
  # Create temporary directory
  temp_dir = tempdir()
  custom_prefix = "test_"
  
  # Check for files with custom names in custom directory
  expected_files = c(
    file.path(temp_dir, paste0(custom_prefix, "varimpByFold.tex")),
    file.path(temp_dir, paste0(custom_prefix, "varimpAll.tex")),
    file.path(temp_dir, paste0(custom_prefix, "varimpConsistent.tex"))
  )
  
  # Ensure cleanup happens even if test fails
  on.exit({
    cleanup_latex_files(dir = temp_dir, outname = custom_prefix, verbose = FALSE)
  })
  
  # Test with custom outname and directory
  exportLatex(vim, outname = custom_prefix, dir = temp_dir)
  
  existing_custom_files = expected_files[file.exists(expected_files)]
  expect_true(length(existing_custom_files) > 0,
              info = "Custom named files should be created in custom directory")
  
  # Test manual cleanup with custom names
  # Files should be cleaned up by on.exit() handler
})
# The tests above run a real varimpact() fit, so they skip under R CMD check
# and never exercise two branches: every variable significant (no hline
# after the p = 0.05 cut-off) and a non-empty consistent table. A hand-built
# result object reaches both without a fit, and writes into a scratch
# directory so it can run anywhere.
test_that("all-significant results and a consistent table are exported", {
  results_all = data.frame(
    Estimate = c(0.30, 0.20),
    "Adj. p-value" = c(0.001, 0.010),
    check.names = FALSE,
    row.names = c("x1", "x2"))
  mock_vim = list(
    results_by_fold = data.frame(x1 = c("1 vs 3", "1 vs 3"), x2 = c("2 vs 3", "1 vs 3"),
                                 row.names = c("fold 1", "fold 2")),
    results_all = results_all,
    results_consistent = results_all)

  out_dir = tempfile("latex")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  result = exportLatex(mock_vim, dir = out_dir)

  expect_s3_class(result$xtables$consistent, "xtable")
  expect_equal(nrow(result$tables$consistent), 2L)
  for (f in c("varimpByFold.tex", "varimpAll.tex", "varimpConsistent.tex")) {
    expect_true(file.exists(file.path(out_dir, f)), info = f)
  }
  expect_true(any(grepl("consisRes", readLines(file.path(out_dir, "varimpConsistent.tex")))))

  # cleanup_latex_files() reports what it removed, and then that there is
  # nothing left to remove.
  expect_output(removed <- cleanup_latex_files(dir = out_dir, verbose = TRUE),
                "Successfully removed 3 of 3")
  expect_equal(removed, rep(TRUE, 3))
  expect_output(again <- cleanup_latex_files(dir = out_dir, verbose = TRUE),
                "No LaTeX files found")
  expect_identical(again, logical(0))
})
