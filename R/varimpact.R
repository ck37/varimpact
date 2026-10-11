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
#' @param impute Type of missing value imputation to conduct. One of:
#'   "median" (default), "knn", "zero". Note: knn results in the covariate
#'   data being centered and scaled.
#' @param miss.cut eliminates explanatory (X) variables with proportion
#' of missing obs > cut.off
#' @param bins_numeric Number of quantile bins used to discretize a numeric
#'   variable when it serves as the variable of interest, provided it has more
#'   than that many distinct values. Note that within each cross-validation fold
#'   these bins are then further aggregated by penalized histogram density
#'   estimation, in order to avoid small cell sizes (Hubbard, Kennedy & van der
#'   Laan 2018, section 9.8). The number of bins actually compared is therefore
#'   data-driven, and is often substantially smaller than \code{bins_numeric} -
#'   in practice frequently just two.
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
#' @param id Optional cluster identifier(s) for cluster-robust standard
#'   errors, for data in which observations are correlated within groups
#'   (siblings in a family, patients in a hospital, students in a school).
#'   The simplest form is the name(s) of columns in \code{data}, e.g.
#'   \code{id = "family"}, or \code{id = c("school", "neighborhood")} for
#'   multiway clustering (Cameron, Gelbach & Miller 2011); those columns are
#'   used as clusters and are not analyzed as predictors. Alternatively pass
#'   the ids directly: a vector with one value per observation, or a list or
#'   data frame of such vectors. When given: every cluster is kept within a
#'   single CV fold; SuperLearner's internal cross-validation keeps clusters
#'   together; influence-curve variances sum within clusters, across all
#'   folds at once, which follows from the cross-validated TMLE's influence
#'   curve (Zheng & van der Laan 2010; Hubbard, Kherad-Pajouh & van der Laan
#'   2016) with clusters as the independent units, with a small-sample (CR3)
#'   correction for clusters of unequal size (MacKinnon, Nielsen & Webb
#'   2023); and confidence intervals
#'   and p-values use a t distribution, as normal critical values over-reject
#'   when clusters are few (Cameron, Gelbach & Miller 2008; Cameron & Miller
#'   2015). Its degrees of freedom, shown in the \code{df} column of the
#'   results, count the clusters behind the scarcer of the variable's low and
#'   high levels (fold by fold), minus one, and are at most (fewest clusters -
#'   1), following Cameron, Gelbach & Miller (2011) for several clustering
#'   variables. A contrast carried by few clusters has few effective degrees
#'   of freedom (Bell & McCaffrey 2002; Imbens & Kolesar 2016); this count is
#'   a simple, conservative heuristic for them rather than a published rule.
#'   A variable
#'   nested in another (family within site) is covered by clustering on the
#'   outer one and is dropped with a message. The asymptotics rely on the
#'   number of clusters, not observations, so results with few clusters
#'   (fewer than about 40) should be read with caution. Learners that tune
#'   themselves by internal cross-validation, such as \code{SL.glmnet}, still
#'   choose their own folds without regard to clusters; that does not affect
#'   the validity of the results, since the CV-TMLE folds keep clusters
#'   intact, but it can make their tuning optimistic. NULL (the default)
#'   assumes independent observations.
#' @param folds_by With several crossed cluster variables, the name of the
#'   column of \code{id} to form CV folds on. By default folds keep every
#'   cluster variable intact when that is possible, and otherwise fall back to
#'   the coarsest variable with at least \code{V} clusters, with a warning.
#'   With crossed variables, training and validation folds share clusters of
#'   the variables not folded on, so a fold's choice of low and high levels
#'   carries over to its validation data through those clusters, which biases
#'   variables that vary with them (a factor measured at that level, say). The
#'   clustered variance cannot correct this, so name here the variable that
#'   the variables of interest vary with.
#' @param adjustment_exclusions Named list specifying adjustment variables to
#'   exclude when estimating the importance of a given variable. Each name is a
#'   column of \code{data} whose importance is being estimated, and each element
#'   is a character vector of other columns of \code{data} to leave out of its
#'   adjustment set, e.g. \code{list(V1 = c("V2", "V3"))}. Use it to keep a
#'   mediator or collider out of a particular adjustment set. Naming a factor
#'   excludes all of its indicator columns, and a variable's missingness
#'   indicator is excluded along with the variable. Defaults to an empty list,
#'   which adjusts for every other variable as usual.
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
#' Bell, R. M., & McCaffrey, D. F. (2002). \emph{Bias reduction in standard
#' errors for linear regression with multi-stage samples}. Survey Methodology,
#' 28(2), 169-181.
#'
#' Benjamini, Y., & Hochberg, Y. (1995). \emph{Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing}. Journal of the
#' royal statistical society. Series B (Methodological), 289-300.
#'
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2008). \emph{Bootstrap-based
#' improvements for inference with clustered errors}. The Review of Economics
#' and Statistics, 90(3), 414-427.
#'
#' Cameron, A. C., Gelbach, J. B., & Miller, D. L. (2011). \emph{Robust
#' inference with multiway clustering}. Journal of Business & Economic
#' Statistics, 29(2), 238-249.
#'
#' Cameron, A. C., & Miller, D. L. (2015). \emph{A practitioner's guide to
#' cluster-robust inference}. Journal of Human Resources, 50(2), 317-372.
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
#' Imbens, G. W., & Kolesar, M. (2016). \emph{Robust standard errors in small
#' samples: some practical advice}. The Review of Economics and Statistics,
#' 98(4), 701-712.
#'
#' MacKinnon, J. G., Nielsen, M. O., & Webb, M. D. (2023). \emph{Cluster-robust
#' inference: a guide to empirical practice}. Journal of Econometrics, 232(2),
#' 272-299.
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
#' Zheng, W., & van der Laan, M. J. (2010). \emph{Asymptotic theory for
#' cross-validated targeted maximum likelihood estimation}. U.C. Berkeley
#' Division of Biostatistics Working Paper Series, Working Paper 273.
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
#' cleanup_latex_files()
#'
#' # Impute by median rather than knn.
#' \dontrun{
#' vim <- varimpact(Y = Y, data = X[, 1:3], impute = "median")
#' }
#'
#' ####################################
#' # Clustered observations: cluster-robust standard errors.
#' \dontrun{
#' # Suppose the 100 observations are 25 families of four.
#' family <- rep(1:25, each = 4)
#' vim <- varimpact(Y = Y, data = X[, 1:3], id = family)
#'
#' # Or name a column of data, which is then not analyzed as a predictor.
#' vim <- varimpact(Y = Y, data = cbind(X[, 1:3], family), id = "family")
#'
#' # Several clustering variables give multiway clustering, e.g.
#' # varimpact(Y = Y, data = data, id = c("school", "neighborhood"))
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
#' future::plan("multisession", workers = 2)
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
           digits = 4L,
           adjustment_exclusions = list(),
           id = NULL,
           folds_by = NULL) {

  # Time the full function execution.
  time_start = proc.time()

  ######################
  # Argument checks.

  # Handle vector input by converting to data frame
  if (is.vector(data) && !is.list(data)) {
    if (verbose) cat("Converting vector input to data frame.\n")
    data <- data.frame(X1 = data)
  }

  # An id written as data$column or data[["column"]] that evaluates to NULL is
  # almost always a typo such as id = data$famly (a missing column gives NULL),
  # which would otherwise silently turn clustering off. Other NULLs, such as a
  # wrapper passing on its own id = NULL default, just mean no clustering.
  id_expr = substitute(id)
  if (is.null(id) && is.call(id_expr) &&
      as.character(id_expr[[1]])[1] %in% c("$", "[[")) {
    warning("id was supplied but is NULL (for example data$column with a ",
            "misspelled column name), so standard errors are NOT clustered.",
            call. = FALSE)
  }

  # id can name columns of data, e.g. id = "family" or id = c("family",
  # "site"). Those columns become the cluster ids and are removed from data so
  # they are not analyzed as predictors. A character vector as long as the
  # data is taken as the ids themselves, not as column names.
  if (is.character(id) && length(id) < nrow(data)) {
    missing_cols = setdiff(id, colnames(data))
    if (length(missing_cols) > 0L) {
      stop("id names column(s) not found in data: ",
           paste(missing_cols, collapse = ", "), ".")
    }
    id_cols = unique(id)
    id = data[, id_cols, drop = FALSE]
    data = data[, setdiff(colnames(data), id_cols), drop = FALSE]
    if (ncol(data) == 0L) {
      stop("After removing the id column(s), data has no variables left.")
    }
    if (verbose) {
      cat("Clustering on column(s)", paste(id_cols, collapse = ", "),
          "of data; they are not analyzed as predictors.\n")
    }
  }

  # Handle single column case with warning
  if (ncol(data) == 1L) {
    warning("Using single variable for variable importance analysis. Results may be limited.")
    if (verbose) cat("Single variable detected in data.\n")
  }

  # Check the exclusion list once here rather than inside the per-fold loops,
  # so that a typo is reported a single time instead of once per fold.
  check_adjustment_exclusions(adjustment_exclusions, colnames(data))

  # Ensure that Y is numeric; e.g. can't be a factor.
  stopifnot(class(Y) %in% c("numeric", "integer"))

  if (length(Y) != nrow(data)) {
    stop("Y has ", length(Y), " values but data has ", nrow(data), " rows.")
  }

  # Validate the cluster ids once, up front. NULL means no clustering, which
  # includes an id in which every observation is its own cluster.
  cluster_id = check_cluster_id(id, length(Y), verbose = TRUE)
  df_cluster = cluster_df(cluster_id)
  if (!is.null(cluster_id) && verbose) {
    cat("Cluster-robust standard errors using:",
        paste0(names(cluster_id), " (",
               vapply(cluster_id, function(x) length(unique(x)), numeric(1)),
               " clusters)", collapse = ", "),
        "\nt critical values with at most", df_cluster,
        "degrees of freedom.\n")
  }

  # The grouping that CV folds keep intact. Predictors whose values each fall
  # within one such group cannot be compared across folds, so they are set
  # aside here rather than dropping out of the results unexplained.
  fold_groups = cluster_fold_groups(cluster_id, V, folds_by)
  data = drop_cluster_level_columns(data, fold_groups)

  if (family == "binomial" &&
      (min(Y, na.rm = TRUE) < 0 || max(Y, na.rm = TRUE) > 1)) {
    stop("With binomial family Y must be bounded by [0, 1]. Specify family=\"gaussian\" otherwise.")
  }

  if (!family %in% c("binomial", "gaussian")) {
    stop('Family must be either "binomial" or "gaussian".')
  }

  # Without this an unrecognized value falls through every branch in
  # process_numerics() and surfaces much later as
  # "sum(is.na(data.numW)) == 0 is not TRUE", which says nothing about the
  # typo that caused it. "mean" is accepted here so that it still reaches its
  # own not-implemented message rather than being reported as a bad value.
  impute = match.arg(impute, c("median", "knn", "zero", "mean"))

  if (parallel && verbose) {
    cat("Future backend set to the following:\n")
    print(future::plan())
  }

  # Save bounds on the full Y variables for later transformation if Y is not binary.
  if (family == "binomial" || length_unique(Y) == 2L) {
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
  # With clustering, whole clusters are assigned to folds.
  folds = create_cv_folds(V, Y, verbose = verbose, id = fold_groups)

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
                adjustment_exclusions = adjustment_exclusions,
                cluster_id = cluster_id,
                fold_groups = fold_groups,
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
                 verbose_reduction = verbose_reduction,
                 adjustment_exclusions = adjustment_exclusions,
                 cluster_id = cluster_id,
                 fold_groups = fold_groups)

  # Combine the separate continuous and factor results.
  results =
    compile_results(numeric_vims$colnames_numeric,
                    factor_vims$colnames_factor,
                    numeric_vims$vim_numeric,
                    factor_vims$vim_factor,
                    V = V,
                    verbose = verbose,
                    df = df_cluster)

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
                   cv_folds = folds,
                   # Cluster ids actually used (after dropping nested and
                   # singleton variables; NULL without clustering) and the
                   # largest degrees of freedom for the t critical values;
                   # each variable's own df is in results_all$df.
                   cluster_id = cluster_id,
                   cluster_df = df_cluster))

  # Set a custom class so that we can override print and summary.
  class(results) = "varimpact"

  invisible(results)
}
