#' @title Variable importance estimation using causal inference (targeted learning)
#'
#' @description \code{varimpact} returns variable importance statistics ordered
#'   by statistical significance using a combination of data-adaptive target
#'   parameter
#'
#' @details
#' The function performs the following functions.
#'  \enumerate{
#'  \item Drops variables missing > miss.cut of time (tuneable).
#'  \item Separate out covariates into factors and continuous (ordered).
#'  \item Drops variables for which their distribution is uneven  - e.g., all 1
#'  value (tuneable) separately for factors and numeric variables (ADD MORE
#'  DETAIL HERE)
#'  \item Makes dummy variable basis for factors, including naming dummies
#'  to be traceable to original factor variables later.
#'  \item Makes new ordered variable of integers mapped to intervals defined by
#'  deciles for the ordered numeric variables (automatically makes) fewer
#'  categories if original variable has < 10 values.
#'  \item Creates associated list of number of unique values and the list of them
#'  for each variable for use in variable importance part.
#'  \item Makes missing covariate basis for both factors and ordered variables
#'  \item For each variable, after assigning it as A, uses optimal histogram
#'  function to combine values using the distribution of A | Y=1 to avoid very
#'  small cell sizes in distribution of Y vs. A (tuneable) (ADD DETAIL)
#'  \item Uses HOPACH* to cluster variables associated confounder/missingness
#'  basis for W, that uses specified minimum number of adjustment variables.
#'  \item Finds min and max estimate of E(Ya) w.r.t. a. after looping through
#'  all values of A* (after processed by histogram)
#'  \item Returns estimate of E(Ya(max)-Ya(min)) with SE using CV-TMLE.
#' }
#' *HOPACH is "Hierarchical Ordered Partitioning and Collapsing Hybrid"
#'
#' @param Y outcome of interest (numeric vector)
#' @param data data frame of predictor variables of interest for
#' which function returns VIM's. (possibly a matrix?)
#' @param A_names Names of the variables for which we want to estimate importance,
#'  a subset of the data argument.
#' @param V Number of cross-validation folds.
#' @param Q.library library used by SuperLearner for model of outcome
#' versus predictors
#' @param g.library library used by SuperLearner for model of
#' predictor variable of interest versus other predictors
#' @param family family ('binomial' or 'gaussian')
#' @param minYs mininum # of obs with event  - if it is < minYs, skip VIM
#' @param minCell is the cut-off for including a category of A in analysis, and
#'   presents the minumum of cells in a 2x2 table of the indicator of that level
#'   versus outcome, separately by training and validation sample.
#' @param adjust_cutoff Maximum number of adjustment variables during TMLE. If
#'   more than this cutoff varimpact will attempt to reduce the dimensions to
#'   that number (using HOPACH hierarchical clustering). Must not be more than
#'   15 due to HOPACH constraints. Set to NULL to disable any dimension
#'   reduction.
#' @param corthres cut-off correlation with explanatory
#' variable for inclusion of an adjustment variables
#' @param impute Type of missing value imputation to conduct. One of: "zero",
#'   "median", "knn" (default). Note: knn results in the covariate data being centered/scaled.
#' @param miss.cut eliminates explanatory (X) variables with proportion
#' of missing obs > cut.off
#' @param bins_numeric Numbers of bins when discretizing numeric variables.
#' @param quantile_probs_factor Quantiles used to check if factors have
#'   sufficient variation.
#' @param quantile_probs_numeric Quantiles used to check if numerics have
#'   sufficient variation.
#' @param parallel Use parallel processing if a backend is registered; enabled
#'   by default.
#' @param verbose Boolean - if TRUE the method will display more detailed
#'   output.
#' @param verbose_tmle Boolean - if TRUE, will display even more detail on the TMLE
#'   estimation process.
#' @param verbose_reduction Boolean - if TRUE, will display more detail during
#'   variable reduction step (clustering).
#' @param digits Number of digits to round the value labels.
#'
#' @return Results object. TODO: add more detail here.
#'
#' @importFrom stats cor model.matrix na.omit pnorm quantile var
#' @importFrom SuperLearner All
#'
#' @seealso
#' \code{\link[varimpact]{exportLatex}}, \code{\link[varimpact]{print.varimpact}}
#'
#' @encoding utf8
#'
#' @section Authors:
#' Alan E. Hubbard and Chris J. Kennedy, University of California, Berkeley
#'
#'
#' @section References:
#' Benjamini, Y., & Hochberg, Y. (1995). \emph{Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing}. Journal of the
#' royal statistical society. Series B (Methodological), 289-300.
#'
#' Gruber, S., & van der Laan, M. J. (2012). \emph{tmle: An R Package for
#' Targeted Maximum Likelihood Estimation}. Journal of Statistical Software,
#' 51(i13).
#'
#' Hubbard, A. E., Kherad-Pajouh, S., & van der Laan, M. J. (2016).
#' \emph{Statistical Inference for Data Adaptive Target Parameters}. The
#' international journal of biostatistics, 12(1), 3-19.
#'
#' Hubbard, A., Munoz, I. D., Decker, A., Holcomb, J. B., Schreiber, M. A.,
#' Bulger, E. M., ... & Rahbar, M. H. (2013). \emph{Time-Dependent Prediction
#' and Evaluation of Variable Importance Using SuperLearning in High Dimensional
#' Clinical Data}. The journal of trauma and acute care surgery, 75(1 0 1), S53.
#'
#' Hubbard, A. E., & van der Laan, M. J. (2016). \emph{Mining with inference:
#' data-adaptive target parameters (pp. 439-452)}. In P. Buhlmann et al. (Ed.),
#' \emph{Handbook of Big Data}. CRC Press, Taylor & Francis Group, LLC: Boca
#' Raton, FL.
#'
#' van der Laan, M. J. (2006). \emph{Statistical inference for variable
#' importance}. The International Journal of Biostatistics, 2(1).
#'
#' van der Laan, M. J., & Pollard, K. S. (2003). \emph{A new algorithm for
#' hybrid hierarchical clustering with visualization and the bootstrap}. Journal
#' of Statistical Planning and Inference, 117(2), 275-303.
#'
#' van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). \emph{Super
#' learner}. Statistical applications in genetics and molecular biology, 6(1).
#'
#' van der Laan, M. J., & Rose, S. (2011). \emph{Targeted learning: causal
#' inference for observational and experimental data}. Springer Science &
#' Business Media.
#'
#' @examples
#' ####################################
#' # Create test dataset.
#' set.seed(1)
#' N <- 100
#' num_normal <- 5
#' X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
#' Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#' # Add some missing data to X so we can test imputation.
#' for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA
#'
#' ####################################
#' # Basic example
#'
#' vim <- varimpact(Y = Y, data = X[, 1:3])
#' vim
#' vim$results_all
#' exportLatex(vim)
#' 
#' # Clean up LaTeX files
#' suppressWarnings({
#'   file.remove(c("varimpByFold.tex", "varimpAll.tex", "varimpConsistent.tex"))
#' })
#'
#' # Impute by median rather than knn.
#' \dontrun{
#' vim <- varimpact(Y = Y, data = X[, 1:3], impute = "median")
#' }
#'
#' ####################################
#' # Multicore parallel example.
#' \dontrun{
#' # Setup multicore parallelization.
#' library(future)
#' plan("multisession", workers = 2)
#'
#' vim <- varimpact(Y = Y, data = X[, 1:3])
#' }
#'
#' ####################################
#' # Cluster parallel example.
#' \dontrun{
#' cl = parallel::makeCluster(2L)
#' plan(cluster, workers = cl)
#' vim <- varimpact(Y = Y, data = X[, 1:3])
#' parallel::stopCluster(cl)
#' }
#'
#' ####################################
#' # mlbench BreastCancer example.
#' \dontrun{
#' data(BreastCancer, package="mlbench")
#' data <- BreastCancer
#'
#' set.seed(1, "L'Ecuyer-CMRG")
#' # Reduce to a dataset of 100 observations to speed up testing.
#  data = data[sample(nrow(data), 100), ]
#
#' # Create a numeric outcome variable.
#' data$Y <- as.numeric(data$Class == "malignant")
#
#' # Use multicore parallelization to speed up processing.
#' future::plan("multiprocess", workers = 2)
#' vim <- varimpact(Y = data$Y, data = subset(data, select=-c(Y, Class, Id)))
#' }
#'
#' @export
varimpact =
  function(Y,
           data,
           A_names = colnames(data),
           V = 2L,
           #Q.library = c("SL.glmnet", "SL.mean"),
           #g.library = c("SL.glmnet", "SL.mean"),
           Q.library = c("SL.glm", "SL.mean"),
           g.library = c("SL.glm", "SL.mean"),
           #g.library = c("SL.stepAIC"),
           family = "binomial",
           minYs = 15L,
           minCell = 0L,
           adjust_cutoff = 10L,
           corthres = 0.8,
           impute = "median",
           miss.cut = 0.5,
           bins_numeric = 10L,
           quantile_probs_factor = c(0.1, 0.9),
           quantile_probs_numeric = quantile_probs_factor,
           verbose = FALSE,
           verbose_tmle = FALSE,
           verbose_reduction = FALSE,
           parallel = TRUE,
           digits = 4L) {

  # Time the full function execution.
  time_start = proc.time()

  ######################
  # Argument checks.

  # Confirm that data has at least two columns.
  if (ncol(data) < 2L) {
    stop("Data argument must have at least two columns.")
  }

  # Ensure that Y is numeric; e.g. can't be a factor.
  stopifnot(class(Y) %in% c("numeric", "integer"))

  if (family == "binomial" &&
      (min(Y, na.rm = TRUE) < 0 || max(Y, na.rm = TRUE) > 1)) {
    stop("With binomial family Y must be bounded by [0, 1]. Specify family=\"gaussian\" otherwise.")
  }

  if (!family %in% c("binomial", "gaussian")) {
    stop('Family must be either "binomial" or "gaussian".')
  }

  if (parallel && verbose) {
    cat("Future backend set to the following:\n")
    print(future::plan())
  }

  # Save bounds on the full Y variables for later transformation if Y is not binary.
  if (family == "binomial" || length(unique(Y)) == 2) {
    #Qbounds = NULL
    Qbounds = c(0, 1)
  } else {
    # This part is duplicated from the TMLE code in tmle_init_stage1.

    # Define Qbounds just for continuous (non-binary) outcomes.
    Qbounds = range(Y, na.rm = TRUE)
    # Extend bounds 10% beyond the observed range.
    # NOTE: if one of the bounds is zero then it won't be extended.
    Qbounds = Qbounds + 0.1 * c(-abs(Qbounds[1]), abs(Qbounds[2]))
  }

  ########
  # Applied to Explanatory (X) data frame
  sna = sapply(data, sum_na)

  n = nrow(data)

  #######
  # Missing proportion by variable.
  mis.prop = sna / n

  #######
  # Cut-off for eliminating variable for proportion of obs missing.
  data = data[, mis.prop < miss.cut, drop = FALSE]

  # TODO: move this stuff into a separate function.
  if (verbose) cat("Removed", sum(mis.prop >= miss.cut), "variables due to high",
                   "missing value proportion.\n")

  # Separate dataframe into factors-only and numerics-only.
  # Also converts characters to factors automatically.
  separated_data = separate_factors_numerics(data)

  factors = process_factors(separated_data$df_factors,
                            quantile_probs_factor = quantile_probs_factor,
                            miss.cut = miss.cut,
                            verbose = verbose)

  # Pre-process numeric/continuous variables.
  numerics =
    process_numerics(separated_data$df_numerics,
                     quantile_probs_numeric = quantile_probs_numeric,
                     miss.cut = miss.cut,
                     bins_numeric = bins_numeric,
                     impute = impute,
                     verbose = verbose)

  cat("Finished pre-processing variables.\n")

  cat("\nProcessing results:\n")
  cat("- Factor variables:", factors$num_factors, "\n")
  cat("- Numeric variables:", numerics$num_numeric, "\n\n")

  # Create cross-validation folds (2 by default).
  folds = create_cv_folds(V, Y, verbose = verbose)

  # VIM for factors.
  factor_vims =
    vim_factors(Y = Y, numerics = numerics, factors = factors,
                V = V, folds = folds,
                A_names = A_names,
                family = family,
                minCell = minCell,
                minYs = minYs,
                Q.library = Q.library,
                g.library = g.library,
                Qbounds = Qbounds,
                corthres = corthres,
                adjust_cutoff = adjust_cutoff,
                verbose = verbose,
                verbose_tmle = verbose_tmle,
                verbose_reduction = verbose_reduction)

  # Repeat for numerics.
  numeric_vims =
    vim_numerics(Y = Y, numerics = numerics, factors = factors,
                 V = V, folds = folds,
                 A_names = A_names,
                 family = family,
                 minCell = minCell,
                 minYs = minYs,
                 Q.library = Q.library,
                 g.library = g.library,
                 Qbounds = Qbounds,
                 corthres = corthres,
                 adjust_cutoff = adjust_cutoff,
                 verbose = verbose,
                 verbose_tmle = verbose_tmle,
                 verbose_reduction = verbose_reduction)

  # Combine the separate continuous and factor results.
  results =
    compile_results(numeric_vims$colnames_numeric,
                    factor_vims$colnames_factor,
                    numeric_vims$vim_numeric,
                    factor_vims$vim_factor,
                    V = V,
                    verbose = verbose)

  # End timing the full execution.
  time_end = proc.time()

  # Final compilation of results.
  results = c(results,
              # Append additional settings to the results object.
              # TODO: make this a sublist?
              list(V = V,
                   g.library = g.library,
                   Q.library = Q.library,
                   minCell = minCell,
                   minYs = minYs,
                   family = family,
                   datafac.dumW  = factors$datafac.dumW,
                   miss.fac = factors$miss.fac,
                   data.numW = numerics$data.numW,
                   numeric_vims = numeric_vims,
                   factor_vims = factor_vims,
                   impute_info = numerics$impute_info,
                   time = time_end - time_start,
                   cv_folds = folds))

  # Set a custom class so that we can override print and summary.
  class(results) = "varimpact"

  invisible(results)
}
