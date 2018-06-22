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
#' @importFrom future future_lapply
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
#' # Setup multicore parallelization.
#' library(future)
#' plan("multiprocess", workers = 2)
#'
#' vim <- varimpact(Y = Y, data = X[, 1:3])
#' vim
#' vim$results_all
#' exportLatex(vim)
#'
#' # Impute by median rather than knn.
#' \dontrun{
#' vim <- varimpact(Y = Y, data = X[, 1:3], impute = "median")
#' }
#'
#' ####################################
#' # Multicore parallel example.
#' \dontrun{

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
           Q.library = c("SL.glmnet", "SL.mean"),
           g.library = c("SL.glmnet", "SL.mean"),
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

  if (verbose) cat("Removed", sum(mis.prop >= miss.cut), "variables due to high",
                   "missing value proportion.\n")


  # Vector of number of unique values by variable.
  # TODO: run in parallel to support very wide/big datasets.
  # TODO: not used.
  num.values = sapply(data, length_unique)

  # Separate dataframe into factors-only and numerics-only.
  # Also converts characters to factors automatically.
  separated_data = separate_factors_numerics(data)

  # Create a dataframe consisting only of columns that are factors.
  data.fac = separated_data$df_factors
  # And samesies for numerics.
  data.num = separated_data$df_numerics


  factors = process_factors(data.fac,
                            quantile_probs_factor = quantile_probs_factor,
                            miss.cut = miss.cut,
                            verbose = verbose)

  ###########################################################
  # Pre-process numeric/continuous variables.

  # TODO: need to convert numeric processing code to its own function so that
  # we can test it formally.
  # preprocess_numerics(data.num, quantile_probs, bins_numeric, verbose = verbose)
  # Return results:
  # which variables were dropped, and why.
  # data.cont.dist -
  # impute info

  numerics =
    process_numerics(data.num,
                     quantile_probs_numeric = quantile_probs_numeric,
                     miss.cut = miss.cut,
                     bins_numeric = bins_numeric,
                     impute = impute,
                     verbose = verbose)

  cat("Finished pre-processing variables.\n")

  # Create cross-validation folds (2 by default).
  folds = create_cv_folds(V, Y, verbose = verbose)

  n = length(Y)

  # TODO: print variable summary before we begin analysis.
  cat("\nProcessing results:\n")
  cat("- Factor variables:", factors$num_factors, "\n")
  cat("- Numeric variables:", numerics$num_numeric, "\n\n")

  ###############################################################################
  # VIM for Factors
  # NOTE: we use && so that conditional will short-circuit if num_factors == 0.
  if (factors$num_factors > 0L && ncol(factors$data.fac) > 0L) {
    cat("Estimating variable importance for", factors$num_factors, "factors.\n")

    # Find the level of covariate that has lowest risk
    datafac.dumW = factors$datafac.dum
    # NOTE: can't we skip this line because we already imputed missing data to 0?
    datafac.dumW[is.na(factors$datafac.dum)] = 0

    #############################
    # Below is to get indexing vectors so that any basis functions related to current A
    # that are in covariate matrix can be removed.
    names.fac = colnames(factors$data.fac)
    nmes.facW = colnames(datafac.dumW)
    nmes.mfacW = colnames(factors$miss.fac)
    nchar.facW = nchar(nmes.facW) + 1
    nchar.mfacW = nchar(nmes.mfacW) + 1

    XXm = regexpr("XX", nmes.facW)
    XXm[XXm < 0] = nchar.facW[XXm < 0]

    XXm2 = regexpr("XX", nmes.mfacW)
    XXm2[XXm2 < 0] = nchar.mfacW[XXm2 < 0]

    vars.facW = substr(nmes.facW, 1, XXm - 1)
    vars.mfacW = substr(nmes.mfacW, 7, XXm2 - 1)

    xc = ncol(data.fac)
    n.fac = nrow(data.fac)

    # vim_factor = lapply(1:xc, function(i) {

    # vim_factor will be a list of results, one element per factor variable.
    # Define var_i just to avoid automated NOTEs, will be overwritten by foreach.
    var_i = NULL
    #vim_factor = foreach::foreach(var_i = 1:xc, .verbose = verbose, .errorhandling = "stop") %do_op% {
    vim_factor = future.apply::future_lapply(1:xc, future.seed = TRUE, function(var_i) {
    #vim_factor = lapply(1:xc, function(var_i) {
      nameA = names.fac[var_i]

      if (verbose) cat("Var:", nameA, var_i, "out of", xc, "factor variables\n")

      if (!nameA %in% A_names) {
        if (verbose) cat("Skipping", nameA, "as it is not in A_names.\n")
        return(NULL)
      }

      # Loop over each fold.
      # TODO: incorporate this loop into parallelization.
      # for (fold_k in 1:V) {

      # This is looping sequentially for now.
      #fold_results = foreach::foreach(fold_k = 1:V) foreach::`%do%` {
      fold_results = lapply(1:V, function(fold_k) {
        if (verbose) cat("i =", var_i, "V =", fold_k, "\n")

        # All data not in this fold is the training data.
        At = data.fac[folds != fold_k, var_i]

        # All data in this fold is the validation data.
        Av = data.fac[folds == fold_k, var_i]

        Yt = Y[folds != fold_k]
        Yv = Y[folds == fold_k]


        #######################################
        # Create adjustment dataframe.

        ### acit.numW is just same as acit.cont.dist except with NA's replaced by
        ### 0's.
        mtch1 = match(vars.facW, nameA)
        mtch2 = match(vars.mfacW, nameA)
        Adum = data.frame(factors$datafac.dum[, is.na(mtch1) == FALSE])
        dumW = factors$datafac.dum[, is.na(mtch1)]
        missdumW = factors$miss.fac[, is.na(mtch2)]

        if (is.null(missdumW)) {
          missdumW = rep(NA, n.fac)
        }
        if (is.null(numerics$miss.cont)) {
          numerics$miss.cont = rep(NA, n.fac)
        }
        if (is.null(dumW)) {
          dumW = rep(NA, n.fac)
        }
        if (is.null(numerics$data.numW)) {
          numerics$data.numW = rep(NA, n.fac)
        }

        W = data.frame(numerics$data.numW, numerics$miss.cont, dumW, missdumW)

        # Restrict to columns in which there is less than 100% missingness.
        W = W[, !apply(is.na(W), 2, all), drop = FALSE]

        #######################################

        # Divide into training and validation subsets.
        Wt = W[folds != fold_k, , drop = FALSE]
        Wv = W[folds == fold_k, , drop = FALSE]

        Adum = data.frame(Adum[folds != fold_k, ])

        ###
        # Pull out any variables that are overly correlated with At (corr coef < corthes)
        #if (sd(Adum) == 0) {
        #  if (verbose) cat("Warning: sd of Adum = 0.\n")
        #}

        # Suppress possible "the standard deviation is zero" warning from cor().
        # TODO: investigate more and confirm that this is ok.
        suppressWarnings({
          corAt = apply(stats::cor(Adum, Wt, use = "complete.obs"), 2, max_sqr)
        })
        corAt[corAt < -1] = 0
        # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
        incc = abs(corAt) < corthres & !is.na(corAt)

        if (verbose && sum(!incc) > 0) {
          cat("Removed", sum(!incc), "columns based on correlation threshold", corthres, "\n")
        }

        Wv = Wv[, incc, drop = F]
        Wt = Wt[, incc, drop = F]

        if (verbose) {
          cat("Columns:", ncol(Wt))
          if (!is.null(adjust_cutoff)) cat(" Reducing dimensions to", adjust_cutoff)
          cat("\n")
        }

        # Use HOPACH to reduce dimension of W to some level of tree
        reduced_results = reduce_dimensions(Wt, Wv, adjust_cutoff, verbose = verbose_reduction)

        Wtsht = reduced_results$data
        Wvsht = reduced_results$newX

        # We should have no constant columns after calling reduce_dimensions().
        # Remove any NA values - but shouldn't these already be imputed?
        is_constant = sapply(Wtsht, function(col) var(col, na.rm = TRUE) == 0)
        # Restrict to just the TRUE variables - those that are constant.
        is_constant = is_constant[is_constant]

        if (verbose) {
          cat("Updated ncols -- training:", ncol(Wtsht), "test:", ncol(Wvsht), "\n")

          # We should have no constant columns after calling reduce_dimensions().
          if (length(is_constant) > 0) {
            cat("Constant columns (", length(is_constant), "):\n")
            print(is_constant)
          }
        }

        # Finished with any needed clustering for variable reduction.

        deltat = as.numeric(!is.na(Yt) & !is.na(At))
        deltav = as.numeric(!is.na(Yv) & !is.na(Av))

        # TODO (CK): don't do this, in order to use the delta missingness estimation.
        # To avoid crashing TMLE function just drop obs missing A or Y if the
        # total number of missing is < 10
        if (sum(deltat == 0) < 10) {
          Yt = Yt[deltat == 1]
          At = At[deltat == 1]
          Wtsht = Wtsht[deltat == 1, , drop = F]
          deltat = deltat[deltat == 1]
        }

        levA = levels(At)

        if (length(unique(Yt)) == 2) {
          # Binary outcome.

          # Minimum numer of observations for each cell in validation fold.
          minc = apply(table(Av, Yv), 1, min)

          # Minimum numer of observations for each cell in training fold.
          minc2 = apply(table(At, Yt), 1, min)

          # Restrict analysis to levels of this variable in which
          # there are sufficient observations in each fold
          # across each treatment and outcome combination.
          vals = levA[pmin(minc, minc2) > minCell]
        } else {
          # Continuous outcome.
          vals = levA
        }
        num.cat = length(vals)

        # CK TODO 6/6: don't assume that positive outcome is the rare
        #  outcome. (e.g. via table)

        # Number of positive outcomes in training data.
        nYt = sum(Yt[!is.na(At)])

        # Number of positive outcomes in validation data.
        nYv = sum(Yv[!is.na(Av)])

        # Create a list to hold the results we calculate in this fold.
        # Set them to default values and update as they are calculated.
        fold_result = list(
          failed = TRUE,
          # Message should report on the status for this fold.
          message = "",
          obs_training = length(Yt),
          obs_validation = length(Yv),
          error_count = 0,
          # Save the variables chosen by reduce_dimensions().
          variables = reduced_results$variables,
          # Results for estimating the maximum level / treatment.
          level_max = list(
            # Level is which bin was chosen.
            level = NULL,
            # Label is the description of that bin.
            label = NULL,
            # val_preds contains the g, Q, and H predictions on the validation data.
            val_preds = NULL,
            # Estimate of EY on the training data.
            estimate_training = NULL,
            # Risk from SuperLearner on Q.
            risk_Q = NULL,
            # Risk from SuperLearner on g.
            risk_g = NULL
          )
        )
        # Copy the blank result to a second element for the minimum level/bin.
        fold_result$level_min = fold_result$level_max

        ############################
        # Don't do if 1) no more than one category of A left or
        # 2) if missingness pattern for A is such that there are few death events left
        # in either (< minYs)
        # Applies only to binary outcomes, not continuous.
        if ((length(unique(Yt)) == 2 && (num.cat < 2 || min(nYt, nYv) < minYs)) ||
            (length(is_constant) > 0 && mean(is_constant) == 1)) {
          if (length(is_constant) > 0 && mean(is_constant) == 1) {
            error_msg = paste("Skipping", nameA, "because HOPACH reduced W to",
                              "all constant columns.")
          } else if (num.cat < 2) {
            error_msg = paste("Skipping", nameA, "due to lack of variation.")
          } else {
            error_msg = paste("Skipping", nameA, "due to minY constraint.", min(nYt, nYv), "<", minYs)
          }
          if (verbose) cat(error_msg, "\n")

          fold_result$message = error_msg
          # At this point here we are skipping to the end of the loop.

        } else {
          if (verbose) cat("Estimating TMLE on training", paste0("(", num.cat, ")"))

          error_count = 0

          training_estimates = list()

          # Estimate Y_a, Q hat and g hat for each level of our current variable,
          # on the training data.
          for (j in 1:num.cat) {

            # Create a treatment indicator, where 1 = obs in this bin
            # and 0 = obs not in this bin.
            IA = as.numeric(At == vals[j])

            # Any observations missing At are assigned to 0.
            IA[is.na(IA)] = 0

            # if(min(table(IA,Yt))>=)

            # CV-TMLE: we are using this for three reasons:
            # 1. Estimate Y_a on training data.
            # 2. Estimate Q on training data.
            # 3. Estimate g on training data.
            tmle_result = try(estimate_tmle2(Yt, IA, Wtsht, family, deltat,
                                    Q.lib = Q.library,
                                    Qbounds = Qbounds,
                                    g.lib = g.library, verbose = verbose_tmle),
                      silent = !verbose)

            if (class(tmle_result) == "try-error") {
              # TMLE estimation failed.
              if (verbose) cat("X")
              error_count = error_count + 1
            } else {
              # TMLE estimation successed.

              # Save label
              tmle_result$label = vals[j]

              training_estimates[[j]] = tmle_result

              if (verbose) {
                cat(".")
              }
            }
          }
          # Finished looping over each level of the assignment variable.
          if (verbose) cat(" done. Errors:", error_count, "\n")

          fold_result$error_count = error_count

          # Extract theta estimates.
          theta_estimates = sapply(training_estimates, function(result) {
            # Handle errors in the tmle estimation by returning NA.
            ifelse("theta" %in% names(result), result$theta, NA)
          })

          if (!all(is.na(theta_estimates))) {
            # Identify maximum EY1 (theta)
            # Note: this may be NA if the tmle estimation failed.
            maxj = which.max(theta_estimates)

            # Identify minimum EY1 (theta)
            # Note: this may be NA if the tmle estimation failed.
            minj = which.min(theta_estimates)
            if (verbose) {
              cat("Max level:", vals[maxj], paste0("(", maxj, ")"),
                    "Min level:", vals[minj], paste0("(", minj, ")"), "\n")
            }
          } else {
            maxj = NA
            minj = NA
          }

          # This fold failed if we got an error for each category
          # Or if the minimum and maximum bin is the same.
          if (error_count == num.cat ||
              (is.na(minj) && is.na(maxj)) ||
              minj == maxj) {
            message = paste("Fold", fold_k, "failed,")
            if (length(theta_estimates) == 0 || error_count == num.cat) {
              message = paste(message, "all", num.cat, "levels had errors.")
            } else {
              message = paste(message, "min and max level are the same. (j = ", minj,
                              "label = ", training_estimates[[minj]]$label, ")")
            }
            fold_result$message = message

            if (verbose) {
              cat(message, "\n")
            }
          } else {

            # Extract max items.
            maxEY1 = training_estimates[[maxj]]$theta
            labmax = vals[maxj]

            # Save these items into the fold_result list.
            fold_result$level_max$level = maxj
            fold_result$level_max$estimate_training = maxEY1
            fold_result$level_max$label = labmax
            # Save the Q risk for the discrete SuperLearner.
            # We don't have the CV.SL results for the full SuperLearner as it's too
            # computationallity intensive.
            fold_result$level_max$risk_Q =
                training_estimates[[maxj]]$q_model$cvRisk[
                  which.min(training_estimates[[maxj]]$q_model$cvRisk)]
            # And the g's discrete SL risk.
            fold_result$level_max$risk_g =
                training_estimates[[maxj]]$g_model$cvRisk[
                  which.min(training_estimates[[maxj]]$g_model$cvRisk)]

            # Extact TMLE results.
            fold_result$level_max

            #fold_result$level_max$tmle = training_estimates[[maxj]]

            # Extract min items.
            minEY1 = training_estimates[[minj]]$theta
            labmin = vals[minj]

            # Save these items into the fold_result list.
            fold_result$level_min$level = minj
            fold_result$level_min$estimate_training = minEY1
            fold_result$level_min$label = labmin
            #fold_result$level_min$tmle = training_estimates[[minj]]

            # Save the Q risk for the discrete SuperLearner.
            # We don't have the CV.SL results for the full SuperLearner as it's too
            # computationallity intensive.
            fold_result$level_min$risk_Q =
              training_estimates[[minj]]$q_model$cvRisk[
                which.min(training_estimates[[minj]]$q_model$cvRisk)]
            # And the g's discrete SL risk.
            fold_result$level_min$risk_g =
              training_estimates[[minj]]$g_model$cvRisk[
                which.min(training_estimates[[minj]]$g_model$cvRisk)]

            # Turn to validation data.

            # Estimate minimum level (control).

            # Indicator for having the desired control bin on validation.
            IA = as.numeric(Av == vals[minj])

            # Missing values are not taken to be in this level.
            IA[is.na(IA)] = 0

            if (verbose) cat("\nMin level prediction - apply_tmle_to_validation()\n")

            # CV-TMLE: predict g, Q, and clever covariate on validation data.
            min_preds = try(apply_tmle_to_validation(Yv, IA, Wvsht, family,
                                     deltav, training_estimates[[minj]],
                                     verbose = verbose))

            # Old version:
            #res = try(estimate_tmle(Yv, IA, Wvsht, family, deltav,
            #                        Q.lib = Q.library,
            #                        g.lib = g.library, verbose = verbose),
            #          silent = T)

            if (class(min_preds) == "try-error") {
              message = paste("CV-TMLE prediction on validation failed during",
                            "low/control level.")
              fold_result$message = message
              if (verbose) cat(message, "\n")
            } else {
              # Save the result.
              fold_result$level_min$val_preds = min_preds

              # Switch to maximum level (treatment).

              # Indicator for having the desired treatment bin on validation
              IA = as.numeric(Av == vals[maxj])

              # Missing values are not taken to be in this level.
              IA[is.na(IA)] = 0

              if (verbose) cat("\nMax level prediction - apply_tmle_to_validation()\n")

              # CV-TMLE: predict g, Q, and clever covariate on validation data.
              max_preds = try(apply_tmle_to_validation(Yv, IA, Wvsht, family,
                                        deltav, training_estimates[[maxj]],
                                        verbose = verbose))
             # Old code:
              #res2 = try(estimate_tmle(Yv, IA, Wvsht, family, deltav,
              #                         Q.lib = Q.library,
              #                         g.lib = g.library, verbose = verbose),
              #           silent = !verbose)


              if (class(max_preds) == "try-error") {
                message = paste("CV-TMLE prediction on validation failed",
                      "during high/treatment level.")
                fold_result$message = message
                if (verbose) cat(message, "\n")
              } else {
                # Save the result.
                fold_result$level_max$val_preds = max_preds
                fold_result$message = "Succcess"
                # Fold succeeded.
                fold_result$failed = FALSE
              }
            }
          }
        }
        if (verbose) cat("Completed fold", fold_k, "\n\n")

        # Return results for this fold.
        fold_result
      }
      ) # End lapply if we're not using foreach.
      # Done looping over each fold.

      # Create list to save results for this variable.
      var_results = list(
        EY1V = NULL,
        EY0V = NULL,
        thetaV = NULL,
        varICV = NULL,
        labV = NULL,
        nV = NULL,
        fold_results = fold_results,
        type = "factor"
      )

      # TODO: compile results into the new estimate.

      if (verbose) cat("Estimating pooled min.\n")
      pooled_min = estimate_pooled_results(lapply(fold_results, function(x) x$level_min),
                                           verbose = verbose)
      cat("\n")
      if (verbose) cat("Estimating pooled max.\n")
      pooled_max = estimate_pooled_results(lapply(fold_results, function(x) x$level_max),
                                           verbose = verbose)
      cat("\n")

      var_results$EY0V = pooled_min$thetas
      var_results$EY1V = pooled_max$thetas

      if (length(var_results$EY1V) == length(var_results$EY0V)) {
        var_results$thetaV = var_results$EY1V - var_results$EY0V
      } else {
        if (verbose) {
          cat("Error: EY1V and EY0V are different lengths. EY1V =",
              length(var_results$EY1V), "EY0V =", length(var_results$EY0V), "\n")
        }
        var_results$thetaV = rep(NA, max(length(var_results$EY1V),
                                         length(var_results$EY0V)))
      }


      # Save how many observations were in each validation fold.
      var_results$nV = sapply(fold_results, function(x) x$obs_validation)

      # Combine labels into a two-column matrix.
      # First column is min and second is max.
      # TODO: not sure if data structure for this part is correct.
      labels = do.call(rbind,
                 lapply(fold_results, function(x) c(x$level_min$label, x$level_max$label)))

      var_results$labV = labels

      # If either of the thetas is null it means that all CV-TMLE folds failed.
      if (!is.null(pooled_min$thetas)) {

        # Influence_curves here is a list, with an element for each fold.
        var_results$varICV = sapply(1:V, function(index) {
          if (length(pooled_max$influence_curves) >= index &&
              length(pooled_min$influence_curves) >= index) {
            var(pooled_max$influence_curves[[index]] - pooled_min$influence_curves[[index]])
          } else {
            NA
          }
        })

        if (verbose) {
          signif_digits = 4

          ey0_mean = mean(pooled_min$thetas)
          if (is.numeric(ey0_mean)) {
            cat("[Min] EY0:", signif(ey0_mean, signif_digits))
            if (is.numeric(pooled_min$epsilon)) {
              cat(" Epsilon:", signif(pooled_min$epsilon, signif_digits))
            }
            cat("\n")
          }

          ey1_mean =  mean(pooled_max$thetas)
          if (is.numeric(ey1_mean)) {
            cat("[Max] EY1:", signif(ey1_mean, signif_digits))
            if (is.numeric(pooled_max$epsilon)) {
              cat(" Epsilon:", signif(pooled_max$epsilon, signif_digits))
            }
            cat("\n")
          }

          cat("ATEs:", signif(var_results$thetaV, signif_digits), "\n")
          cat("Variances:", signif(var_results$varICV, signif_digits), "\n")
          cat("Labels:\n")
          print(labels)
          cat("\n")
        }
      }

      # Return results for this factor variable.
      var_results
    #} # End foreach loop over all variables.
    }) # End lapply or future_lapply if we're not using foreach.

    if (verbose) cat("Factor VIMs:", length(vim_factor), "\n\n")

    # Confirm that we have the correct number of results, otherwise fail out.
    stopifnot(length(vim_factor) == xc)

    colnames_factor = colnames(data.fac)
  } else {
    colnames_factor = NULL
    vim_factor = NULL
    cat("No factor variables - skip VIM estimation.\n\n")
  }

  ###############################################################################
  # Repeat for Continuous Explanatory Variables
  ###############################################################################

  cor.two = function(x, y) {
    (stats::cor(na.omit(cbind(x, y)))[1, 2])^2
  }

  ######################################################
  # We use && so that the second check will be skipped when num_numeric == 0.
  if (numerics$num_numeric > 0 && ncol(numerics$data.cont.dist) > 0) {
    cat("Estimating variable importance for", numerics$num_numeric, "numerics.\n")

    xc = ncol(numerics$data.cont.dist)
    names.cont = colnames(numerics$data.cont.dist)
    n.cont = nrow(numerics$data.cont.dist)

    # Tally the number of unique values (bins) in each numeric variable; save as a vector.
    numcat.cont = apply(numerics$data.cont.dist, 2, length_unique)

    cats.cont = lapply(1:xc, function(i) {
      sort(unique(numerics$data.cont.dist[, i]))
    })

    ### Loop over each numeric variable.
    # Define var_i just to avoid automated NOTEs, will be overwritten by foreach.
    var_i = NULL
    #vim_numeric = foreach::foreach(var_i = 1:num_numeric, .verbose = verbose,
    #                               .errorhandling = "stop") %do_op% {
    vim_numeric = future.apply::future_lapply(1:numerics$num_numeric, future.seed = TRUE,
                                        function(var_i) {
    #vim_numeric = lapply(1:num_numeric, function(var_i) {
      nameA = names.cont[var_i]

      if (verbose) cat("i =", var_i, "Var =", nameA, "out of", xc, "numeric variables\n")

      if (!nameA %in% A_names) {
        if (verbose) cat("Skipping", nameA, " as it is not in A_names.\n")
        return(NULL)
      }

      #for (fold_k in 1:V) {
      # This is looping sequentially for now.
      #fold_results = foreach (fold_k = 1:V) %do% {
      # TODO: convert to future_lapply
      fold_results = lapply(1:V, function(fold_k) {
        if (verbose) cat("Fold", fold_k, "of", V, "\n")

        # data.cont.dist is the discretized version of the numeric variables of size bins_numeric
        At = numerics$data.cont.dist[folds != fold_k, var_i]
        Av = numerics$data.cont.dist[folds == fold_k, var_i]
        Yt = Y[folds != fold_k]
        Yv = Y[folds == fold_k]

        # Create a list to hold the results we calculate in this fold.
        # Set them to default values and update as they are calculated.
        fold_result = list(
          # Set failed = FALSE at the very end if everything works.
          failed = TRUE,
          # Message should report on the status for this fold.
          message = "",
          obs_training = length(Yt),
          obs_validation = length(Yv),
          error_count = 0,
          # Results for estimating the maximum level / treatment.
          level_max = list(
            # Level is which bin was chosen.
            level = NULL,
            # Label is the string description of that bin.
            label = NULL,
            # val_preds contains the g, Q, and H predictions on the validation data.
            val_preds = NULL,
            # Estimate of EY on the training data.
            estimate_training = NULL,
            # Risk from SuperLearner on Q.
            risk_Q = NULL,
            # Risk from SuperLearner on g.
            risk_g = NULL
          )
        )
        # Copy the blank result to a second element for the minimum level/bin.
        fold_result$level_min = fold_result$level_max

        # Conduct penalized histogramming by looking at the distribution of the rare outcome over
        # the treatment variable. So we create A_Y1 as the conditional distribution of treatment given Y = 1.
        # P(A | Y = 1).
        if (length(unique(Yt)) == 2L) {
          # Binary outcome.

          A_Y1 = At[Yt == 1 & !is.na(At)]

          # Check if AY1 has only a single value. If so, skip histogramming to avoid an error.
          singleAY1 = length(unique(na.omit(A_Y1))) == 1L
        } else {
          # Continuous outcome - just restricted to non-missing treatment.
          A_Y1 = At[!is.na(At)]
          singleAY1 = F
        }

        if (!singleAY1) {
          # Within this CV-TMLE fold look at further combining bins of the treatment based on the
          # penalized histogram.
          # Note that this is only examining the already discretized version of the treatment variable.
          penalized_hist = histogram::histogram(A_Y1, verbose = F, type = "irregular", plot = F)
          hh = penalized_hist$breaks

          # TODO: see if these next two steps are ever used/needed.

          # Check if the final cut-off is less that the maximum possible level; if so extend to slightly
          # larger than the maximimum possible level.
          if (hh[length(hh)] < max(At, na.rm = TRUE)) {
            hh[length(hh)] = max(At, na.rm = TRUE) + 0.1
          }

          # Check if the lowest cut-off is greater than the minimum possible bin; if so extend to slightly
          # below the minimum level.
          if (hh[1] > min(At[At > 0], na.rm = TRUE)) {
            hh[1] = min(At[At > 0], na.rm = TRUE) - 0.1
          }

          # Re-bin the training and validation vectors for the treatment variable based on the penalized
          # histogram.
          # This is creating factors, with levels specific to this CV-TMLE fold.
          Atnew = cut(At, breaks = hh)
          Avnew = cut(Av, breaks = hh)

          # TODO: check if the binning results in no-variation, and handle separately from the below situation.

        }
        if (singleAY1 || length(na.omit(unique(Atnew))) <= 1 ||
            length(na.omit(unique(Avnew))) <= 1) {
          error_msg = paste("Skipping", nameA, "in this fold because there is no variation.")
          if (verbose) cat(error_msg, "\n")
          fold_result$message = error_msg
          #warning(error_msg)
        } else {

          # These labels are simply the quantiles right now.
          At_bin_labels = names(table(Atnew))

          # Non-discretized version of A in the training data; converted to a vector.
          At_raw = data.num[folds != fold_k, var_i]

          # Loop over the Atnew levels and figure out the equivalent true range of this bin
          # by examining the non-discretized continuous variable.
          for (newlevel_i in 1:length(unique(as.numeric(na.omit(Atnew))))) {
            range = range(na.omit(At_raw[na.omit(which(as.numeric(Atnew) == newlevel_i))]))
            label_i = paste0("[", round(range[1], 2), ", ", round(range[2], 2), "]")
            At_bin_labels[newlevel_i] = label_i
          }
          At_bin_labels

          Atnew = as.numeric(Atnew) - 1
          Avnew = as.numeric(Avnew) - 1

          # Update the number of bins for this numeric variable.
          # CK: note though, this is specific to this CV-TMLE fold - don't we
          # need to differentiate which fold we're in?
          numcat.cont[var_i] = length(At_bin_labels)

          # change this to match what was done for factors - once
          # cats.cont[[i]]=as.numeric(na.omit(unique(Atnew)))
          cats.cont[[var_i]] = as.numeric(names(table(Atnew)))

          ### acit.numW is just same as data.cont.dist except with NA's replaced by
          ### 0's.
          # TODO: cbind iteratively to create W matrix below, so that we don't
          # need these extra NA vectors.
          if (is.null(numerics$miss.cont)) {
            numerics$miss.cont = rep(NA, n.cont)
          }
          if (is.null(factors$miss.fac)) {
            factors$miss.fac = rep(NA, n.cont)
          }
          if (is.null(factors$datafac.dumW)) {
            factors$datafac.dumW = rep(NA, n.cont)
          }
          if (is.null(numerics$data.numW)) {
            numerics$data.numW = rep(NA, n.cont)
          }

          # Construct a matrix of adjustment variables in which we use the imputed dataset
          # but remove the current treatment variable.
          W = data.frame(numerics$data.numW[, -var_i, drop = FALSE],
                         numerics$miss.cont,
                         factors$datafac.dumW,
                         factors$miss.fac)

          # Remove any columns in which all values are NA.
          # CK: but we're using imputed data, so there should be no NAs actually.
          # (With the exception of the NA vectors possibly added above.
          W = W[, !apply(is.na(W), 2, all), drop = FALSE]

          # Separate adjustment matrix into the training and test folds.
          Wt = W[folds != fold_k, , drop = FALSE]
          Wv = W[folds == fold_k, , drop = FALSE]

          # Identify the missingness indicator for this treatment.
          miss_ind_name = paste0("Imiss_", nameA)

          # Remove the missingness indicator for this treatment (if it exists) from the adjustment set.
          Wt = Wt[, colnames(Wt) != miss_ind_name, drop = FALSE]
          Wv = Wv[, colnames(Wt) != miss_ind_name, drop = FALSE]

          # Pull out any variables that are overly correlated with At (corr coef > corthes)

          # Suppress possible warning from cor() "the standard deviation is zero".
          # TODO: remove those constant variables beforehand?
          suppressWarnings({
            corAt = apply(Wt, 2, cor.two, y = At)
          })


          keep_vars = corAt < corthres & !is.na(corAt)

          if (verbose && sum(!keep_vars) > 0) {
            cat("Removed", sum(!keep_vars), "columns based on correlation threshold", corthres, "\n")
          }

          Wv = Wv[, keep_vars, drop = FALSE]
          Wt = Wt[, keep_vars, drop = FALSE]

          if (verbose) {
            cat("Columns:", ncol(Wt))
            if (!is.null(adjust_cutoff)) cat(" Reducing dimensions to", adjust_cutoff)
            cat("\n")
          }

          # Use HOPACH to reduce dimension of W to some level of tree.
          reduced_results = reduce_dimensions(Wt, Wv, adjust_cutoff, verbose = verbose_reduction)

          Wtsht = reduced_results$data
          Wvsht = reduced_results$newX

          # Identify any constant columns.
          is_constant = sapply(Wtsht, function(col) var(col) == 0)
          is_constant = is_constant[is_constant]

          if (verbose) {
            cat("Updated ncols -- training:", ncol(Wtsht), "test:", ncol(Wvsht), "\n")
            # Restrict to true elements.
            if (length(is_constant) > 0) {
              cat("Constant columns (", length(is_constant), "):\n")
              print(is_constant)
            }
          }

          # Indicator that Y and A are both defined.
          deltat = as.numeric(!is.na(Yt) & !is.na(Atnew))
          deltav = as.numeric(!is.na(Yv) & !is.na(Avnew))

          # TODO: may want to remove this procedure, which is pretty arbitrary.
          if (sum(deltat == 0) < 10) {
            Yt = Yt[deltat == 1]
            Wtsht = Wtsht[deltat == 1, , drop = F]
            Atnew = Atnew[deltat == 1]
            deltat = deltat[deltat == 1]
          }

          vals = cats.cont[[var_i]]

          Atnew[is.na(Atnew)] = -1
          Avnew[is.na(Avnew)] = -1

          if ((length(is_constant) > 0 && mean(is_constant) == 1) ||
              (length(unique(Yt)) == 2 && min(table(Avnew[Avnew >= 0], Yv[Avnew >= 0])) <= minCell)) {
            if (length(is_constant) > 0 && mean(is_constant) == 1) {
              error_msg = paste("Skipping", nameA, "because HOPACH reduced W",
                "to all constant columns.")
            } else {
              error_msg = paste("Skipping", nameA, "due to minCell constraint.\n")
            }
            if (T || verbose) cat(error_msg)
            fold_result$message = error_msg
            # warning(error_msg)
            # Go to the next loop iteration.
            #next
          } else {
            # CK TODO: this is not exactly the opposite of the IF above. Is that intentional?
            #if (length(unique(Yt)) > 2 || min(table(Avnew, Yv)) > minCell) {

            # Tally how many bins fail with an error.
            error_count = 0

            if (verbose) cat("Estimating training TMLEs", paste0("(", numcat.cont[var_i], " bins)"))

            training_estimates = list()

            # Loop over each bin for this variable.
            for (j in 1:numcat.cont[var_i]) {

              # Create a treatment indicator, where 1 = obs in this bin
              # and 0 = obs not in this bin.
              IA = as.numeric(Atnew == vals[j])

              # CV-TMLE: we are using this for three reasons:
              # 1. Estimate Y_a on training data.
              # 2. Estimate Q on training data.
              # 3. Estimate g on training data.
              tmle_result = try(estimate_tmle2(Yt, IA, Wtsht, family, deltat,
                                               Q.lib = Q.library,
                                               # Pass in Q bounds from the full
                                               # range of Y (training & test).
                                               Qbounds = Qbounds,
                                               g.lib = g.library, verbose = verbose_tmle),
                                silent = !verbose)

              # Old way:
              #res = try(estimate_tmle(Yt, IA, Wtsht, family, deltat, Q.lib = Q.library,
              #                        g.lib = g.library, verbose = verbose), silent = T)

              if (class(tmle_result) == "try-error") {
                # Error.
                if (verbose) cat("X")
                error_count = error_count + 1
              } else {

                # TMLE succeeded (hopefully).

                # Save bin label.
                tmle_result$label = At_bin_labels[j]

                training_estimates[[j]] = tmle_result

                if (verbose) cat(".")
              }
            }
            # Finished looping over each level of the assignment variable.
            if (verbose) cat(" done.\n")

            fold_result$error_count = error_count

            # Extract theta estimates.
            theta_estimates = sapply(training_estimates, function(result) {
              # Handle errors in the tmle estimation by returning NA.
              ifelse("theta" %in% names(result), result$theta, NA)
            })

            # Identify maximum EY1 (theta)
            maxj = which.max(theta_estimates)

            # Identify minimum EY1 (theta)
            minj = which.min(theta_estimates)

            if (verbose) {
              cat("Max level:", vals[maxj], At_bin_labels[maxj], paste0("(", maxj, ")"),
                  "Min level:", vals[minj], At_bin_labels[minj], paste0("(", minj, ")"), "\n")
            }

            # Save that estimate.
            maxEY1 = training_estimates[[maxj]]$theta
            labmax = vals[maxj]

            # Save these items into the fold_result list.
            fold_result$level_max$level = maxj
            fold_result$level_max$estimate_training = maxEY1
            #fold_result$level_max$label = labmax
            fold_result$level_max$label = At_bin_labels[maxj]

            # Save the Q risk for the discrete SuperLearner.
            # We don't have the CV.SL results for the full SuperLearner as it's too
            # computationallity intensive.
            fold_result$level_max$risk_Q =
              training_estimates[[maxj]]$q_model$cvRisk[
                which.min(training_estimates[[maxj]]$q_model$cvRisk)]
            # And the g's discrete SL risk.
            fold_result$level_max$risk_g =
              training_estimates[[maxj]]$g_model$cvRisk[
                which.min(training_estimates[[maxj]]$g_model$cvRisk)]


            minEY1 = training_estimates[[minj]]$theta
            labmin = vals[minj]

            # Save these items into the fold_result list.
            fold_result$level_min$level = minj
            fold_result$level_min$estimate_training = minEY1
            #fold_result$level_min$label = labmin
            fold_result$level_min$label = At_bin_labels[minj]

            # Save the Q risk for the discrete SuperLearner.
            # We don't have the CV.SL results for the full SuperLearner as it's too
            # computationallity intensive.
            fold_result$level_min$risk_Q =
              training_estimates[[minj]]$q_model$cvRisk[
                which.min(training_estimates[[minj]]$q_model$cvRisk)]
            # And the g's discrete SL risk.
            fold_result$level_min$risk_g =
              training_estimates[[minj]]$g_model$cvRisk[
                which.min(training_estimates[[minj]]$g_model$cvRisk)]

            # This fold failed if we got an error for each category
            # Or if the minimum and maximum bin is the same.
            if (error_count == numcat.cont[var_i] || minj == maxj) {
              message = paste("Fold", fold_k, "failed,")
              if (error_count == numcat.cont[var_i]) {
                message = paste(message, "all", num.cat, "levels had errors.")
              } else {
                message = paste(message, "min and max level are the same. (j = ", minj,
                                "label = ", training_estimates[[minj]]$label, ")")
              }
              fold_result$message = message

              if (verbose) {
                cat(message, "\n")
              }
            } else {

              # Turn to validation data.

              # Estimate minimum level (control).

              # Indicator for having the desired control bin on validation.
              IA = as.numeric(Avnew == vals[minj])

              # Missing values are not taken to be in this level.
              IA[is.na(IA)] = 0

              if (verbose) cat("\nMin level prediction - apply_tmle_to_validation()\n")

              # CV-TMLE: predict g, Q, and clever covariate on validation data.
              min_preds = try(apply_tmle_to_validation(Yv, IA, Wvsht, family,
                                                       deltav, training_estimates[[minj]],
                                                       verbose = verbose))

              # Old version:
              #res = try(estimate_tmle(Yv, IA, Wvsht, family, deltav,
              #                        Q.lib = Q.library,
              #                        g.lib = g.library, verbose = verbose),
              #          silent = T)

              if (class(min_preds) == "try-error") {
                message = paste("CV-TMLE prediction on validation failed during",
                                "low/control level.")
                fold_result$message = message
                if (verbose) cat(message, "\n")
              } else {
                # Save the result.
                fold_result$level_min$val_preds = min_preds

                # Switch to maximum level (treatment).

                # Indicator for having the desired treatment bin on validation
                IA = as.numeric(Avnew == vals[maxj])

                # Missing values are not taken to be in this level.
                IA[is.na(IA)] = 0

                if (verbose) cat("\nMax level prediction - apply_tmle_to_validation()\n")

                # CV-TMLE: predict g, Q, and clever covariate on validation data.

                max_preds = try(apply_tmle_to_validation(Yv, IA, Wvsht, family,
                                                         deltav, training_estimates[[maxj]],
                                                         verbose = verbose))
                # Old code:
                #res2 = try(estimate_tmle(Yv, IA, Wvsht, family, deltav,
                #                         Q.lib = Q.library,
                #                         g.lib = g.library, verbose = verbose),
                #           silent = !verbose)


                if (class(max_preds) == "try-error") {
                  message = paste("CV-TMLE prediction on validation failed",
                                  "during high/treatment level.")
                  fold_result$message = message
                  if (verbose) cat(message, "\n")
                } else {
                  # Save the result.
                  fold_result$level_max$val_preds = max_preds
                  fold_result$message = "Succcess"
                  fold_result$failed = FALSE
                }
              }
            }
          }
        }
        if (verbose) cat("Completed fold", fold_k, "\n\n")

        # Return results for this fold.
        fold_result
      }) # End lapply
      # Done looping over each fold.

      # Create list to save results for this variable.
      var_results = list(
        EY1V = NULL,
        EY0V = NULL,
        thetaV = NULL,
        varICV = NULL,
        labV = NULL,
        nV = NULL,
        fold_results = fold_results,
        type = "factor",
        name = nameA
      )


      # TODO: compile results into the new estimate.

      # if (verbose) cat("Estimating pooled min.\n")
      pooled_min = estimate_pooled_results(lapply(fold_results, function(x) x$level_min),
                                           verbose = verbose)
      # if (verbose) cat("Estimating pooled max.\n")
      pooled_max = estimate_pooled_results(lapply(fold_results, function(x) x$level_max),
                                           verbose = verbose)

      var_results$EY0V = pooled_min$thetas
      var_results$EY1V = pooled_max$thetas
      if (length(var_results$EY1V) == length(var_results$EY0V)) {
        var_results$thetaV = var_results$EY1V - var_results$EY0V
      } else {
        if (verbose) {
          cat("Error: EY1V and EY0V are different lengths. EY1V =",
              length(var_results$EY1V), "EY0V =", length(var_results$EY0V), "\n")
        }
        var_results$thetaV = rep(NA, max(length(var_results$EY1V),
                                         length(var_results$EY0V)))
      }


      # Save how many observations were in each validation fold.
      var_results$nV = sapply(fold_results, function(x) x$obs_validation)

      # Combine labels into a two-column matrix.
      # First column is min and second is max.
      # TODO: not sure if data structure for this part is correct.
      labels = do.call(rbind,
                       lapply(fold_results, function(x) c(x$level_min$label, x$level_max$label)))

      var_results$labV = labels

      # If either of the thetas is null it means that all CV-TMLE folds failed.
      if (!is.null(pooled_min$thetas)) {

        # Influence_curves here is a list, with each element a set of results.
        var_results$varICV = sapply(1:V, function(index) {
          if (length(pooled_max$influence_curves) >= index &&
              length(pooled_min$influence_curves) >= index) {
            var(pooled_max$influence_curves[[index]] - pooled_min$influence_curves[[index]])
          } else {
            NA
          }
        })


        if (verbose) {
          signif_digits = 4
          ey0_mean = mean(pooled_min$thetas)
          if (is.numeric(ey0_mean)) {
            cat("[Min] EY0:", signif(ey0_mean, signif_digits))
            if (is.numeric(pooled_min$epsilon)) {
              cat(" Epsilon:", signif(pooled_min$epsilon, signif_digits), "\n")
            }
          }

          ey1_mean = mean(pooled_max$thetas)
          if (is.numeric(ey1_mean)) {
            cat("[Max] EY1:", signif(ey1_mean, signif_digits))
            if (is.numeric(pooled_max$epsilon)) {
              cat(" Epsilon:", signif(pooled_max$epsilon, signif_digits), "\n")
            }
          }

          cat("ATEs:", signif(var_results$thetaV, signif_digits), "\n")
          cat("Variances:", signif(var_results$varICV, signif_digits), "\n")
          cat("Labels:\n")
          print(labels)
          cat("\n")
        }

      }

      # Return results for this factor variable.
      var_results
    #} # end foreach loop.
    }) # End lapply or future_lapply if we're not using foreach

    if (verbose) cat("Numeric VIMs:", length(vim_numeric), "\n")

    # Confirm that we have the correct number of results, otherwise fail out.
    if (length(vim_numeric) != numerics$num_numeric) {
      # TODO: remove this.
      save(vim_numeric, file = "varimpact.RData")
      # TEMP remove this:
      stop(paste("We have", numerics$num_numeric, "continuous variables but only",
                 length(vim_numeric), "results."))
    }

    colnames_numeric = colnames(numerics$data.cont.dist)
  } else {
    colnames_numeric = NULL
    vim_numeric = NULL
    data.numW = NULL
    cat("No numeric variables for variable importance estimation.\n")
  }

  if (verbose) cat("Completed numeric variable importance estimation.\n")


  #####################################################
  # Combine the separate continuous and factor results.

  results = compile_results(colnames_numeric,
                            colnames_factor,
                            vim_numeric,
                            vim_factor,
                            V = V,
                            verbose = verbose)

  # End timing the full execution.
  time_end = proc.time()

  # Final compilation of results.
  results = c(results,
              # Append additional settings to the results object.
              # TODO: make this a sublist?
              list(V = V, g.library = g.library, Q.library = Q.library,
                   minCell = minCell, minYs = minYs,
                   family = family,  datafac.dumW  = datafac.dumW,
                   miss.fac = factors$miss.fac,
                   data.numW = data.numW, impute_info = numerics$impute_info,
                   time = time_end - time_start,
                   cv_folds = folds))


  # Set a custom class so that we can override print and summary.
  class(results) = "varimpact"

  invisible(results)
}
