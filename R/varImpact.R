#' @title Variable importance estimation using causal inference (TMLE)
#'
#' @description \code{varImpact} returns variable importance statistics ordered
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
#'  \item Removes spaces from factor labels (used for naming dummies later).
#'  \item Removes spaces from variable names.
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
#'   more than this cutoff varImpact will attempt to reduce the dimensions to
#'   that number (using HOPACH hierarchical clustering). Must not be more than
#'   15 due to HOPACH constraints. Set to NULL to disable any dimension
#'   reduction.
#' @param corthres cut-off correlation with explanatory
#' variable for inclusion of an adjustment variables
#' @param impute Type of missing value imputation to conduct. One of: "zero",
#'   "median", "knn" (default).
#' @param miss.cut eliminates explanatory (X) variables with proportion
#' of missing obs > cut.off
#' @param parallel Use parallel processing if a backend is registered; enabled
#'   by default.
#' @param verbose Boolean - if TRUE the method will display more detailed
#'   output.
#'@param verbose_tmle Boolean - if TRUE, will display even more detail on the TMLE
#'   estimation process.
#' @param digits Number of digits to round the value labels.
#'
#' @return Results object.
#'
#' @importFrom stats cor model.matrix na.omit pnorm quantile var
#'
#' @seealso
#' \code{\link[varImpact]{exportLatex}}, \code{\link[varImpact]{print.varImpact}}
#'
#' @section Authors:
#' Alan E. Hubbard and Chris J. Kennedy, University of California, Berkeley
#'
#' @encoding utf8
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
#' Hubbard, A., Munoz, I. D., Decker, A., Holcomb, J. B., Schreiber, M. A.,
#' Bulger, E. M., ... & Rahbar, M. H. (2013). \emph{Time-Dependent Prediction
#' and Evaluation of Variable Importance Using SuperLearning in High Dimensional
#' Clinical Data}. The journal of trauma and acute care surgery, 75(1 0 1), S53.
#'
#' Hubbard, A. E., & van der Laan, M. J. (2016). \emph{Mining with inference:
#' data-adaptive target parameters (pp. 439-452)}. In P. Bühlmann et al. (Ed.),
#' \emph{Handbook of Big Data}. CRC Press, Taylor & Francis Group, LLC: Boca
#' Raton, FL.
#'
#' van der Laan, M. J. (2006). Statistical inference for variable importance.
#' The International Journal of Biostatistics, 2(1).
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
#' num_normal <- 7
#' X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
#' Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
#' # Add some missing data to X so we can test imputation.
#' for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA
#'
#' ####################################
#' # Basic example
#' vim <- varImpact(Y = Y, data = X)
#' vim
#' vim$results_all
#' exportLatex(vim)
#'
#' # Impute by median rather than knn.
#' vim <- varImpact(Y = Y, data = X, impute = "median")
#'
#' ####################################
#' # doMC parallel (multicore) example.
#' \dontrun{
#' library(doMC)
#' registerDoMC()
#' # Use L'Ecuyer for multicore seeds; see ?set.seed for details.
#' set.seed(23432, "L'Ecuyer-CMRG")
#' vim <- varImpact(Y = Y, data = X)
#' }
#'
#' ####################################
#' # doSNOW parallel example.
#' \dontrun{
#' library(doSNOW)
#' library(RhpcBLASctl)
#' # Detect the number of physical cores on this computer using RhpcBLASctl.
#' cluster <- makeCluster(get_num_cores())
#' registerDoSNOW(cluster)
#' vim <- varImpact(Y = Y, data = X)
#' stopCluster(cluster)
#' }
#'
#' ####################################
#' # mlbench BreastCancer example.
#' data(BreastCancer, package="mlbench")
#' data <- BreastCancer
#
#' # Create a numeric outcome variable.
#' data$Y <- as.numeric(data$Class == "malignant")
#
#' # Use multicore parallelization to speed up processing.
#' \dontrun{
#' doMC::registerDoMC()
#' }
#' vim <- varImpact(Y = data$Y, data = subset(data, select=-c(Y, Class, Id)))
#'
#' @export
varImpact = function(Y,
                     data,
                     A_names = colnames(data),
                     V = 2,
                     Q.library = c("SL.glmnet", "SL.mean"),
                     g.library = c("SL.stepAIC"),
                     family = "binomial",
                     minYs = 15,
                     minCell = 0,
                     adjust_cutoff = 10,
                     corthres = 0.8,
                     impute = "knn",
                     miss.cut = 0.5,
                     verbose = F,
                     verbose_tmle = F,
                     parallel = T,
                     digits = 4) {
  # We need to explictly load SuperLearner due to an issue with
  # the "All" screener as of 2016-12-06.
  # CRAN checks want us to use requireNamespace() instead of library()
  # requireNamespace("SuperLearner")

  # Time the full function execution.
  time_start = proc.time()

  # Ensure that Y is numeric; e.g. can't be a factor.
  stopifnot(class(Y) %in% c("numeric", "integer"))

  if (family == "binomial" && (min(Y, na.rm = T) < 0 || max(Y, na.rm = T) > 1)) {
    stop("With binomial family Y must be bounded by [0, 1]. Specify family=\"gaussian\" otherwise.")
  }

  if (!family %in% c("binomial", "gaussian")) {
    stop('Family must be either "binomial" or "gaussian".')
  }

  # Setup parallelism. Thanks to Jeremy Coyle's origami package for this approach.
  `%do_op%` = foreach::`%do%`
  # Use parallelism if there is a backend registered, unless parallel == F.
  if (foreach::getDoParRegistered() && parallel) {
    `%do_op%` = foreach::`%dopar%`
    if (verbose) cat("Parallel backend detected: using foreach parallelization.\n")
  } else {
    if (verbose) cat("No parallel backend detected. Operating sequentially.\n")
  }


  ########
  # Applied to Explanatory (X) data frame
  sna = apply(data, 2, sum_na)

  n = nrow(data)

  #######
  # Missing proportion by variable.
  mis.prop = sna / n

  #######
  # Cut-off for eliminating variable for proportion of obs missing.
  data = data[, mis.prop < miss.cut]

  if (verbose) cat("Removed", sum(mis.prop >= miss.cut), "variables due to high",
                   "missing value proportion.\n")

  # Function that counts # of unique values.
  length_unique = function(x) {
    length(unique(x))
  }

  # Vector of number of unique values by variable.
  # TODO: run in parallel to support very wide/big datasets.
  num.values = apply(data, 2, length_unique)

  # Separate dataframe into factors-only and numerics-only.
  # Also converts characters to factors automatically.
  separated_data = separate_factors_numerics(data)

  # Create a dataframe consisting only of columns that are factors.
  data.fac = separated_data$df_factors
  # And samesies for numerics.
  data.num = separated_data$df_numerics

  #####################
  if (ncol(data.fac) > 0) {

    if (verbose) cat("Processing factors. Start count:", ncol(data.fac), "\n")

    ######################################
    # Replace blank factor values with NA's.

    # We re-use this num_cols variable in the next section.
    num_cols = ncol(data.fac)
    for (i in 1:num_cols) {
      new_factor = as.character(data.fac[, i])
      # The exclude argument replaces any empty strings with NAs.
      new_factor = factor(new_factor, exclude = "")
      data.fac[, i] = new_factor
    }

    ###################
    # For each factor, apply function and get rid of those where
    # 'true' data.fac is data frame of variables that are factors

    data.fac = restrict_by_quantiles(data.fac, quantile_probs = c(0.1, 0.9))

    dropped_cols = num_cols - ncol(data.fac)

    if (verbose) {
      if (dropped_cols > 0) {
        cat("Dropped", dropped_cols, "factors due to lack of variation.\n")
      } else {
        cat("No factors dropped due to lack of variation.\n")
      }
    }

    # We don't seem to use this yet.
    num.cat = apply(data.fac, 2, length_unique)

    ######################
    # Remove columns with missing data % greater than the threshold.
    sum_nas = apply(data.fac, 2, sum_na)

    if (verbose) cat("Factors with missingness:", sum(sum_nas > 0), "\n")

    miss_pct = sum_nas / nrow(data.fac)

    data.fac = data.fac[, miss_pct < miss.cut, drop = F]

    if (verbose) {
      cat("Dropped", sum(miss_pct >= miss.cut), "factors due to the missingness threshold.\n")
    }

    # Save how many separate factors we have in this dataframe.
    num_factors = ncol(data.fac)

    factor_results = factors_to_indicators(data.fac, verbose = verbose)

    datafac.dum = factor_results$data
    # Here 1 = defined, 0 = missing.
    miss.fac = factor_results$missing_indicators

    if (verbose) {
      cat("End factor count:", num_factors, "Indicators:", ncol(datafac.dum),
          "Missing indicators:", ncol(miss.fac), "\n")
    }
  } else {
    n.fac = 0
    num_factors = 0
    datafac.dumW = NULL
    miss.fac = NULL
    datafac.dum = NULL
  }

  # Finished pre-processing factor variables.
  ###########################################################

  ###########################################################
  # Pre-process numeric/continuous variables.

  if (ncol(data.num) > 0) {
    num_cols = ncol(data.num)
    if (verbose) cat("Processing numerics. Start count:", num_cols, "\n")

    # Remove columns where the 0.1 and 0.9 quantiles have the same value, i.e. insufficent variation.
    # TODO: set this is a configurable setting?
    data.num = restrict_by_quantiles(data.num, quantile_probs = c(0.1, 0.9))

    if (verbose) {
      num_dropped = num_cols - ncol(data.num)
      if (num_dropped > 0) {
        cat("Dropped", num_dropped, "numerics due to lack of variation.\n")
      } else {
        cat("No numerics dropped due to lack of variation.\n")
      }
    }

    # Save how many numeric variables we have in this dataframe.
    num_numeric = ncol(data.num)
  } else {
    num_numeric = 0
  }

  if (num_numeric > 0) {
    if (verbose) cat("Cleaning up", num_numeric, "numeric variables.\n")
    # Make deciles for continuous variables
    X = data.num
    xc = dim(X)[2]
    qt = apply(na.omit(X), 2, quantile, probs = seq(0.1, 0.9, 0.1))
    newX = NULL
    coln = NULL
    varn = colnames(X)

    num.cat = apply(X, 2, length_unique)

    Xnew = NULL

    for (k in 1:num_numeric) {
      Xt = X[, k]
      # cat("Numeric", k, "", colnames(X)[k], "mean missing:", mean(is.na(Xt)), "\n")

      # Suppress the warning that can occur when there are fewer than 10 bins.
      # We should be able to see this as var_binned containing fewer than 10 columns.
      # Warning is in .cut2(): min(xx[xx > upper])
      # "no non-missing arguments to min; returning Inf"
      suppressWarnings({
        # Discretize into up to 10 deciles.
        # TODO: number of bins should be a function argument.
        var_binned = as.numeric(arules::discretize(Xt,
                                                   method = "frequency",
                                                   categories = 10,
                                                   ordered = T))
      })
      Xnew = cbind(Xnew, var_binned)
    }
    colnames(Xnew) = varn
    data.cont.dist = Xnew

    ###############
    # Missing Basis for numeric variables, post-binning.
    xp = ncol(data.cont.dist)
    n.cont = nrow(data.cont.dist)

    sum_nas = apply(data.cont.dist, 2, sum_na)
    nmesX = colnames(data.cont.dist)
    miss.cont = NULL
    nmesm = NULL

    # Create imputed version of the numeric dataframe.
    data.numW = data.num

    # Loop over each binned numeric variable.
    for (k in 1:xp) {
      # Check if that variable has any missing values.
      if (sum_nas[k] > 0) {
        # The effect is that the basis is set to 1 if it exists and 0 if it's missing.
        ix = as.numeric(!is.na(data.cont.dist[, k]))
        miss.cont = cbind(miss.cont, ix)
        nmesm = c(nmesm, paste("Imiss_", nmesX[k], sep = ""))
      }
    }
    # if(is.null(miss.cont)){miss.cont= rep(1,n.cont)}
    colnames(miss.cont) = nmesm

    # Impute missing data in numeric columns.
    if (impute == "zero") {
      data.numW[is.na(data.num)] = 0
      impute_info = 0
    } else if (impute == "median") {
      impute_info = caret::preProcess(data.num, method = "medianImpute")
      data.numW = predict(impute_info, data.num)
    } else if (impute == "mean") {
      stop("Mean imputation not implemented yet. Please use another imputation method.")
    } else if (impute == "knn") {
      impute_info = caret::preProcess(data.num, method = "knnImpute")
      data.numW = predict(impute_info, data.num)
    }

    # Confirm that there are no missing values remaining in data.numW
    stopifnot(sum(is.na(data.numW)) == 0)
  } else {
    impute_info = NULL
    miss.cont = NULL
  }

  cat("Finished pre-processing variables.\n")

  # Create cross-validation folds (2 by default).
  folds = create_cv_folds(V, Y, verbose = verbose)

  n = length(Y)

  ###############################################################################
  # VIM for Factors
  # NOTE: we use && so that conditional will short-circuit if num_factors == 0.
  if (num_factors > 0 && ncol(data.fac) > 0) {
    cat("Estimating variable importance for", num_factors, "factors.\n")

    # Find the level of covariate that has lowest risk
    datafac.dumW = datafac.dum
    # NOTE: can't we skip this line because we already imputed missing data to 0?
    datafac.dumW[is.na(datafac.dum)] = 0

    #############################
    # Below is to get indexing vectors so that any basis functions related to current A
    # that are in covariate matrix can be removed.
    names.fac = colnames(data.fac)
    nmes.facW = colnames(datafac.dumW)
    nmes.mfacW = colnames(miss.fac)
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
    vim_factor = foreach::foreach(var_i = 1:xc, .verbose = F, .errorhandling = "stop") %do_op% {
    #vim_factor = lapply(1:xc, function(var_i) {
      nameA = names.fac[var_i]

      if (verbose) cat("i =", var_i, "Var =", nameA, "out of", xc, "factor variables\n")

      if (!nameA %in% A_names) {
        if (verbose) cat("Skipping", nameA, " as it is not in A_names.\n")
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

        ### acit.numW is just same as acit.cont.dist except with NA's replaced by
        ### 0's.
        mtch1 = match(vars.facW, nameA)
        mtch2 = match(vars.mfacW, nameA)
        Adum = data.frame(datafac.dum[, is.na(mtch1) == F])
        dumW = datafac.dum[, is.na(mtch1)]
        missdumW = miss.fac[, is.na(mtch2)]

        if (!exists("missdumW") || is.null(missdumW)) {
          missdumW = rep(NA, n.fac)
        }
        if (!exists("miss.cont") || is.null(miss.cont)) {
          miss.cont = rep(NA, n.fac)
        }
        if (!exists("dumW") || is.null(dumW)) {
          dumW = rep(NA, n.fac)
        }
        if (!exists("data.numW") || is.null(data.numW)) {
          data.numW = rep(NA, n.fac)
        }

        W = data.frame(data.numW, miss.cont, dumW, missdumW)

        # Restrict to columns in which there is less than 100% missingness.
        W = W[, !apply(is.na(W), 2, all), drop = F]

        # Divide into training and validation subsets.
        Wt = W[folds != fold_k, , drop = F]
        Wv = W[folds == fold_k, , drop = F]

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

        if (verbose) cat("Columns:", ncol(Wt), "Reducing dimensions to", adjust_cutoff, "\n")

        # Use HOPACH to reduce dimension of W to some level of tree
        reduced_results = reduce_dimensions(Wt, Wv, adjust_cutoff, verbose = F)

        Wtsht = reduced_results$data
        Wvsht = reduced_results$newX

        is_constant = sapply(Wtsht, function(col) var(col) == 0)
        is_constant = is_constant[is_constant]

        if (verbose) {
          cat("Updated ncols, training:", ncol(Wtsht), "test:", ncol(Wvsht), "\n")
          # Restrict to true elements.
          if (length(is_constant) > 0) {
            cat("Constant columns (", length(is_constant), "):\n")
            print(is_constant)
          }
        }

        # Finished with any needed clustering for variable reduction.

        deltat = as.numeric(!is.na(Yt) & !is.na(At))
        deltav = as.numeric(!is.na(Yv) & !is.na(Av))

        # To avoid crashing TMLE function just drop obs missing A or Y if the
        # total number of missing is < 10
        if (sum(deltat == 0) < 10) {
          Yt = Yt[deltat == 1]
          At = At[deltat == 1]
          Wtsht = Wtsht[deltat == 1, ]
          deltat = deltat[deltat == 1]
        }

        levA = levels(At)

        if (length(unique(Yt)) == 2) {
          # Binary outcome.

          # Minimum numer of observations in validation fold.
          minc = apply(table(Av, Yv), 1, min)

          # Minimum numer of observations in training fold.
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
          # Message should report on the status for this fold.
          message = "",
          obs_training = length(Yt),
          obs_validation = length(Yv),
          error_count = 0,
          # Results for estimating the maximum level / treatment.
          level_max = list(
            # Level is which bin was chosen.
            level = NULL,
            # Label is the description of that bin.
            label = NULL,
            # val_preds contains the g, Q, and H predictions on the validation data.
            val_preds = NULL,
            # Estimate of EY on the training data.
            estimate_training = NULL
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
                                    Qbounds = range(Y),
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
            if (verbose) cat("maxj:", maxj, "minj:", minj, "\n")
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
            #fold_result$level_max$tmle = training_estimates[[maxj]]

            # Extract min items.
            minEY1 = training_estimates[[minj]]$theta
            labmin = vals[minj]

            # Save these items into the fold_result list.
            fold_result$level_min$level = minj
            fold_result$level_min$estimate_training = minEY1
            fold_result$level_min$label = labmin
            #fold_result$level_min$tmle = training_estimates[[minj]]


            # Turn to validation data.

            # Estimate minimum level (control).

            # Indicator for having the desired control bin on validation.
            IA = as.numeric(Av == vals[minj])

            # Missing values are not taken to be in this level.
            IA[is.na(IA)] = 0

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
              }
            }
          }
        }
        if (verbose) cat("Completed fold", fold_k, "\n")
        fold_result$message = "Succcess"

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
        }
      }

      # Return results for this factor variable.
      var_results
    } # End foreach loop over all variables.
    #) # End lapply if we're not using foreach. (temporary tweak)

    if (verbose) cat("Factor VIMs:", length(vim_factor), "\n")

    # Confirm that we have the correct number of results, otherwise fail out.
    stopifnot(length(vim_factor) == xc)

    colnames_factor = colnames(data.fac)
  } else {
    colnames_factor = NULL
    vim_factor = NULL
    cat("No factor variables - skip VIM estimation.\n")
  }

  ###############################################################################
  # Repeat for Continuous Explanatory Variables
  ###############################################################################

  cor.two = function(x, y) {
    (stats::cor(na.omit(cbind(x, y)))[1, 2])^2
  }

  ######################################################
  # We use && so that the second check will be skipped when num_numeric == 0.
  if (num_numeric > 0 && ncol(data.cont.dist) > 0) {
    cat("Estimating variable importance for", num_numeric, "numerics.\n")

    xc = ncol(data.cont.dist)
    names.cont = colnames(data.cont.dist)
    n.cont = nrow(data.cont.dist)

    # TODO: describe this stpe.
    numcat.cont = apply(data.cont.dist, 2, length_unique)
    # CK: I think we can comment out this line:
    #xc = length(numcat.cont)
    cat("xc is", xc, "and length numcat.cont is", length(numcat.cont), "\n")

    cats.cont = lapply(1:xc, function(i) {
      sort(unique(data.cont.dist[, i]))
    })

    ### Loop over each numeric variable.
    #vim_numeric = lapply(1:xc, function(i) {
    vim_numeric = foreach::foreach(var_i = 1:num_numeric) %do_op% {
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
      fold_results = lapply(1:V, function(fold_k) {
        if (verbose) cat("v =", fold_k, "out of V =", V, "\n")

        At = data.cont.dist[folds != fold_k, var_i]
        Av = data.cont.dist[folds == fold_k, var_i]
        Yt = Y[folds != fold_k]
        Yv = Y[folds == fold_k]

        # Create a list to hold the results we calculate in this fold.
        # Set them to default values and update as they are calculated.
        fold_result = list(
          # Message should report on the status for this fold.
          message = "",
          obs_training = length(Yt),
          obs_validation = length(Yv),
          error_count = 0,
          # Results for estimating the maximum level / treatment.
          level_max = list(
            # Level is which bin was chosen.
            level = NULL,
            # Label is the description of that bin.
            label = NULL,
            # val_preds contains the g, Q, and H predictions on the validation data.
            val_preds = NULL,
            # Estimate of EY on the training data.
            estimate_training = NULL
          )
        )
        # Copy the blank result to a second element for the minimum level/bin.
        fold_result$level_min = fold_result$level_max

        if (length(unique(Yt)) == 2) {
          # Binary outcome.
          AY1 = At[Yt == 1 & !is.na(At)]
          # Check if AY1 has only a single value. If so, skip histogramming to avoid an error.
          singleAY1 = length(unique(na.omit(AY1))) == 1
        } else {
          # Continuous outcome.
          AY1 = At[!is.na(At)]
          singleAY1 = F
        }

        if (!singleAY1) {
          hh = histogram::histogram(AY1, verbose = F, type = "irregular", plot = F)$breaks
          if (hh[length(hh)] < max(At, na.rm = T)) {
            hh[length(hh)] = max(At, na.rm = T) + 0.1
          }
          if (hh[1] > min(At[At > 0], na.rm = T)) {
            hh[1] = min(At[At > 0], na.rm = T) - 0.1
          }
          Atnew = cut(At, breaks = hh)
          Avnew = cut(Av, breaks = hh)
        }
        if (singleAY1 || length(na.omit(unique(Atnew))) <= 1 ||
            length(na.omit(unique(Avnew))) <= 1) {
          error_msg = paste("Skipping", nameA, "in this fold because there is no variation.")
          if (verbose) cat(error_msg, "\n")
          fold_result$message = error_msg
          #warning(error_msg)
        } else {
          #if (length(na.omit(unique(Atnew))) > 1 & length(na.omit(unique(Avnew))) > 1) {
          labs = names(table(Atnew))
          Atnew = as.numeric(Atnew) - 1
          Avnew = as.numeric(Avnew) - 1
          numcat.cont[var_i] = length(labs)
          # change this to match what was done for factors - once
          # cats.cont[[i]]=as.numeric(na.omit(unique(Atnew)))
          cats.cont[[var_i]] = as.numeric(names(table(Atnew)))
          ### acit.numW is just same as data.cont.dist except with NA's replaced by
          ### 0's.
          if (!exists("miss.cont") || is.null(miss.cont)) {
            miss.cont = rep(NA, n.cont)
          }
          if (!exists("miss.fac") || is.null(miss.fac)) {
            miss.fac = rep(NA, n.cont)
          }
          if (!exists("datafac.dumW") || is.null(datafac.dumW)) {
            datafac.dumW = rep(NA, n.cont)
          }
          if (!exists("data.numW") || is.null(data.numW)) {
            data.numW = rep(NA, n.cont)
          }
          W = data.frame(data.numW[, -var_i, drop = F], miss.cont, datafac.dumW, miss.fac)
          W = W[, !apply(is.na(W), 2, all), drop = F]
          Wt = W[folds != fold_k, , drop = F]
          Wv = W[folds == fold_k, , drop = F]

          nmesW = names(Wt)
          mtch = match(nmesW, paste0("Imiss_", nameA))
          Wt = Wt[, is.na(mtch), drop = F]
          Wv = Wv[, is.na(mtch), drop = F]
          ### Pull out any variables that are overly correlated with At (corr coef
          ### < corthes)
          # Suppress possible warning from cor() "the standard deviation is zero".
          suppressWarnings({
            corAt = apply(Wt, 2, cor.two, y = At)
          })
          # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
          incc = corAt < corthres & is.na(corAt) == F
          Wv = Wv[, incc, drop = F]
          Wt = Wt[, incc, drop = F]

          if (verbose) cat("Columns:", ncol(Wt), "Reducing dimensions to", adjust_cutoff, "\n")

          # Use HOPACH to reduce dimension of W to some level of tree
          reduced_results = reduce_dimensions(Wt, Wv, adjust_cutoff, verbose = F)

          Wtsht = reduced_results$data
          Wvsht = reduced_results$newX

          is_constant = sapply(Wtsht, function(col) var(col) == 0)
          is_constant = is_constant[is_constant]

          if (verbose) {
            cat("Updated ncols, training:", ncol(Wtsht), "test:", ncol(Wvsht), "\n")
            # Restrict to true elements.
            if (length(is_constant) > 0) {
              cat("Constant columns (", length(is_constant), "):\n")
              print(is_constant)
            }
          }

          deltat = as.numeric(is.na(Yt) == F & is.na(Atnew) == F)
          deltav = as.numeric(is.na(Yv) == F & is.na(Avnew) == F)

          if (sum(deltat == 0) < 10) {
            Yt = Yt[deltat == 1]
            Wtsht = Wtsht[deltat == 1, ]
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
                                               Qbounds = range(Y),
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

                # Save bin
                tmle_result$label = labs[j]

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
            # Save that estimate.
            maxEY1 = training_estimates[[maxj]]$theta
            labmax = vals[maxj]

            # Save these items into the fold_result list.
            fold_result$level_max$level = maxj
            fold_result$level_max$estimate_training = maxEY1
            fold_result$level_max$label = labmax

            # Identify minimum EY1 (theta)
            minj = which.min(theta_estimates)
            minEY1 = training_estimates[[minj]]$theta
            labmin = vals[minj]

            # Save these items into the fold_result list.
            fold_result$level_min$level = minj
            fold_result$level_min$estimate_training = minEY1
            fold_result$level_min$label = labmin

            # This fold failed if we got an error for each category
            # Or if the minimum and maximum bin is the same.
            if (error_count == numcat.cont[var_i] || minj == maxj) {
              message = paste("Fold", fold_k, "failed,")
              if (error_count == num.cat) {
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
                }
              }
            }
          }
        }
        cat("Completed fold", fold_k, "\n")
        fold_result$message = "Succcess"

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
        }

      }

      # Return results for this factor variable.
      var_results
    } # end foreach loop.
    #) # End lapply if we're not using foreach

    if (verbose) cat("Numeric VIMs:", length(vim_numeric), "\n")

    # Confirm that we have the correct number of results, otherwise fail out.
    if (length(vim_numeric) != num_numeric) {
      save(vim_numeric, file="varimpact.RData")
      # TEMP remove this:
      stop(paste("We have", num_numeric, "continuous variables but only",
                 length(vim_numeric), "results."))
    }

    colnames_numeric = colnames(data.cont.dist)
  } else {
    colnames_numeric = NULL
    vim_numeric = NULL
    data.numW = NULL
    cat("No numeric variables for variable importance estimation.\n")
  }

  if (verbose) cat("Completed numeric variable importance estimation.\n")


  #####################################################
  # Combine the separate continuous and factor results.

  num_numeric = length(colnames_numeric)
  num_factor = length(colnames_factor)

  variable_types = c(rep("ordered", num_numeric), rep("factor", num_factor))
  variable_names = c(colnames_numeric, colnames_factor)

  vim_combined = c(vim_numeric, vim_factor)
  names(vim_combined) = variable_names

  element_length = sapply(vim_combined, length)

  # May increase this to 8, hence >= operator in following lines.
  expected_length = 7

  # Restrict to vim results that have at least 7 elements.
  vim_combined = vim_combined[element_length >= expected_length]
  variable_names = variable_names[element_length >= expected_length]
  variable_types = variable_types[element_length >= expected_length]

  # Set defaults for variables we want to return.
  outres = outres.all = outres.cons = outres.byV = NULL

  # Get rid of any variables that have a validation sample with no
  # estimates of variable importance.
  if (length(vim_combined) == 0) {
    error_msg = "No variable importance estimates could be calculated due to sample size, etc."
    if (verbose) {
      cat(error_msg, "\n")
      cat("Lengths:", element_length, "\n")
    }
    warning(error_msg)
    # TODO: also write output to the file in a separate function call.
    #write("No VIM's could be calculated due to sample size, etc",
    #     file = "AllReslts.tex")
  } else {
    length_ey1 = sapply(vim_combined, function(x) {
      # In some cases class(x$EY1) is NULL, so this avoids a warning.
      if (is.null(x$EY1V)) {
        return(0)
      }
      ey1_vector = na.omit(x$EY1V)
      length(ey1_vector)
    })

    valid_results = vim_combined[length_ey1 == V]
    variable_names = variable_names[length_ey1 == V]
    variable_types = variable_types[length_ey1 == V]


    if (length(valid_results) == 0) {
      error_msg = "No VIMs could be calculated due to sample size, etc."
      if (verbose) {
        cat(error_msg, "\n")
        cat("EY1 lengths:", length_ey1, "\n")
      }
      warning(error_msg)
    } else {

      # Order of results:
      # 1. EY1V, #2. EY0V, #3. thetaV, #4. varICV, #5, labV
      # #6. nV, #7. type, #8. name.

      theta = do.call(rbind, lapply(valid_results, function(x) x$thetaV))

      theta_sum_na = apply(theta, 1, sum_na)
      results_no_na = valid_results[theta_sum_na == 0]
      variable_names = variable_names[theta_sum_na == 0]
      variable_types = variable_types[theta_sum_na == 0]

      # Extact the various components of the results.
      EY1 = do.call(rbind, lapply(results_no_na, function(x) x$EY1))
      EY0 = do.call(rbind, lapply(results_no_na, function(x) x$EY0))
      theta = do.call(rbind, lapply(results_no_na, function(x) x$thetaV))
      varIC = do.call(rbind, lapply(results_no_na, function(x) x$varICV))
      nV = do.call(rbind, lapply(results_no_na, function(x) x$nV))
      n = sum(nV[1, ])

      SEV = sqrt(varIC / nV)

      # Get labels for each of the training sample
      extract_labels = function(x, total_folds) {
        labels = rep(1:total_folds, 2)
        oo = order(labels)
        labels = labels[oo]
        out = as.vector(t(x))
        names(out) = paste0("v.", labels, rep(c("a_L", "a_H"), total_folds))
        out
      }

      # labV is result element 5.
      tst = lapply(results_no_na, function(x) x$labV)
      tst = lapply(tst, extract_labels, total_folds = V)
      labels = do.call(rbind, tst)

      meanvarIC = apply(varIC, 1, mean)

      psi = apply(theta, 1, mean)
      SE = sqrt(meanvarIC / n)

      ci_lower = psi - 1.96 * SE
      ci_upper = psi + 1.96 * SE

      # Number of significant digits.
      signif_digits = 3

      CI95 = paste0("(", signif(ci_lower, signif_digits), " - ", signif(ci_upper, signif_digits), ")")

      # 1-sided p-value
      pvalue = 1 - pnorm(psi / SE)

      ##### FOR THETA (generalize to chi-square test?)
      # TT = (theta[,1] - theta[,2]) / sqrt(SEV[,1]^2 + SEV[,2]^2)
      # pval.comp=2*(1-pnorm(abs(TT))) FOR levels
      # (just make sure in same order)

      num_continuous = sum(variable_types == "ordered")
      num_vars = length(variable_types)

      length.uniq = function(x) {
        length(unique(x))
      }

      ##################
      # Ordered variables first
      cons = NULL
      if (num_continuous > 0) {
        dir = NULL
        for (i in 1:V) {
          lower_temp = labels[1:num_continuous, i * 2 - 1]
          xx = regexpr(",", lower_temp)
          lwr = as.numeric(substr(lower_temp, 2, xx - 1))

          upper_temp = labels[1:num_continuous, i * 2]
          xx = regexpr(",", upper_temp)
          nx = nchar(upper_temp)
          uwr = as.numeric(substr(upper_temp, xx + 1, nx - 1))

          dir = cbind(dir, uwr > lwr)
        }

        cons = apply(dir, 1, length.uniq)
      }

      ##################
      # Factors
      num_factors = num_vars - num_continuous
      if (num_factors > 0) {
        lwr = NULL
        uwr = NULL
        for (i in 1:V) {
          lwr = cbind(lwr, labels[(num_continuous + 1):num_vars, i * 2 - 1])
          uwr = cbind(uwr, labels[(num_continuous + 1):num_vars, i * 2])
        }
        # Count have many unique levels are used for the lower bin - want it to be 1.
        conslwr = apply(lwr, 1, length.uniq)
        # Count how many unique levels are used for the upper bin - want it to be 1.
        consupr = apply(uwr, 1, length.uniq)
        # If conslwr * consupr == 1 then the variable is consistent.
        cons = c(cons, conslwr * consupr)
      }
      # consist= (cons==1 & pval.comp > 0.05)
      signsum = function(x) {
        sum(sign(x))
      }

      # Consistent results need to have all positive all negative thetas,
      # And use the same factor levels for the low and high bins in each CV fold.
      # CK: but really, shouldn't they all be positive? May want to remove abs()
      consist = cons == 1 & abs(apply(theta, 1, signsum)) == V

      procedures = c("Holm", "BH")
      if (num_vars > 1) {
        # Adjust p-values for multiple testing.
        res = multtest::mt.rawp2adjp(pvalue, procedures)
        sorted_rows = res$index
        # This indexing sorts the results in ascending order of unadjusted p-value.
        outres = data.frame(var_type = variable_types[sorted_rows],
                            theta[sorted_rows, , drop = F],
                            psi[sorted_rows],
                            CI95[sorted_rows],
                            res$adj,
                            labels[sorted_rows, , drop = F],
                            consist[sorted_rows])
      } else if (num_vars == 1) {
        # No need for multiple testing adjustment.
        # TODO: just integrate into previous part?
        outres = data.frame(var_type = variable_types,
                            theta,
                            psi,
                            CI95,
                            rawp = pvalue,
                            Holm = pvalue,
                            BH = pvalue,
                            lbs,
                            consist)
      } else {
        outres = NULL
      }

      # TODO: this will give an error if we have no results.

      # Restrict to variables that aren't missing their p-value.
      outres = outres[!is.na(outres[, "rawp"]), , drop = F]

      names(outres)[1:(1 + 2 * V)] = c("VarType", paste0("psiV", 1:V), "AvePsi", "CI95")
      names(outres)[(9 + 2 * V)] = "Consistent"

      ################
      # Get Consistency Measure and only significant
      # TODO: Make BH cut-off flexible in future versions (default at 0.05)
      outres.cons = outres[outres[, "BH"] < 0.05, , drop = F]
      outres.cons = outres.cons[outres.cons[, "Consistent"],
                                c("VarType", "AvePsi", "rawp", "BH", "CI95"), drop = F]
      colnames(outres.cons) = c("Type", "Estimate", "P-value", "Adj. P-value", "CI 95")


      # drops = c('VarType','description','Holm,')
      # outres.all=outres[,!(names(outres) %in% drops)]
      outres.byV = outres[, c(2:(2 + V - 1), 9:(9 + 2 * V)), drop = F]
      outres.all = outres[, c("VarType", "AvePsi", "CI95", "rawp", "BH", "Consistent"), drop = F]
      colnames(outres.all) = c("Type", "Estimate", "CI95", "P-value", "Adj. p-value", "Consistent")
    }
  }

  # End timing the full execution.
  time_end = proc.time()

  # Return results.
  results = list(results_consistent = outres.cons,
                 results_all = outres.all,
                 results_by_fold = outres.byV,
                 results_raw = outres,
                 V = V, g.library = g.library, Q.library = Q.library,
                 minCell = minCell, minYs = minYs,
                 family = family,  datafac.dumW  = datafac.dumW,
                 miss.fac = miss.fac,
                 data.numW = data.numW, impute_info = impute_info,
                 time = time_end - time_start)
  # Set a custom class so that we can override print and summary.
  class(results) = "varImpact"
  invisible(results)
}
