#' @title Variable Impact Estimation
#'
#' @description \code{varImpact} returns variable importance statistics ordered
#' by statistical significance using a combination of data-adaptive target parameter
#'
#' @details
#' The function performs the following functions.
#'  \enumerate{
#'  \item Drops variables missing > miss.cut of time (tuneable).
#'  \item Separate out covariates into factors and continuous (ordered).
#'  \item Drops variables for which their distribution is uneven  - e.g., all 1 value (tuneable)
#'  separately for factors and numeric variables (ADD MORE DETAIL HERE)
#'  \item Changes all factors to remove spaces (used for naming dummies later)
#'  \item Changes variable names to remove spaces
#'  \item Makes dummy variable basis for factors, including naming dummies
#'  to be traceable to original factor variable laters
#'  \item Makes new ordered variable of integers mapped to intervals defined by deciles for the ordered numeric variables (automatically makes)
#'  fewer categories if original variable has < 10 values.
#'  \item Creates associated list of number of unique values and the list of them
#'  for each variable for use in variable importance part.
#'  \item Makes missing covariate basis for both factors and ordered variables
#'  \item For each variable, after assigning it as A, uses
#'  optimal histogram function to combine values using the
#'  distribution of A | Y=1 to avoid very small cell sizes in
#'  distribution of Y vs. A (tuneable) (ADD DETAIL)
#'  \item Uses HOPACH to cluster variables associated confounder/missingness basis for W,
#'  that uses specified minimum number of adjustment variables.
#'  \item Finds min and max estimate of E(Ya) w.r.t. a. after looping through
#'  all values of A* (after processed by histogram)
#'  \item Returns estimate of E(Ya(max)-Ya(min)) with SE
#'  \item Things to do include implementing CV-TMLE and allow reporting of results
#'  that randomly do not have estimates for some of validation samples.
#' }
#'
#' @param Y outcome of interest (numeric vector)
#' @param data data frame of predictor variables of interest for
#' which function returns VIM's. (possibly a matrix?)
#' @param V Number of cross-validation folds.
#' @param Q.library library used by SuperLearner for model of outcome
#' versus predictors
#' @param g.library library used by SuperLearner for model of
#' predictor variable of interest versus other predictors
#' @param family family ('binomial' or 'gaussian')
#' @param minYs mininum # of obs with event  - if it is < minYs, skip VIM
#' @param minCell is the cut-off for including a category of
#' A in analysis, and  presents the minumum of cells in a 2x2 table of the indicator of
#' that level versus outcome, separately by training and validation
#' sample
#' @param ncov minimum number of covariates to include as adjustment variables (must
#' be less than # of basis functions of adjustment matrix)
#' @param corthres cut-off correlation with explanatory
#' variable for inclusion of an adjustment variables
#' @param impute Type of missing value imputation to conduct. One of: "zero",
#'   "median", "knn" (default).
#' @param miss.cut eliminates explanatory (X) variables with proportion
#' of missing obs > cut.off
#' @param parallel Use parallel processing if a backend is registered; enabled by default.
#' @param verbose Boolean - if TRUE the method will display more detailed output.
#'
#' @return Results object.
#'
#' @importFrom stats cor model.matrix na.omit pnorm quantile var
#'
#' @seealso
#' \code{\link[varImpact]{exportLatex}}, \code{\link[varImpact]{print.varImpact}} method
#'
#' @section Authors:
#' Alan E. Hubbard and Chris J. Kennedy, University of California, Berkeley
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
#' Hubbard, A. E., & van der Laan, M. J. (2016). \emph{Mining with inference:
#' data-adaptive target parameter (pp. 439-452)}. In P. Bühlmann et al. (Ed.),
#' \emph{Handbook of Big Data}. CRC Press, Taylor & Francis Group, LLC: Boca
#' Raton, FL.
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
#' N <- 200
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
#' library(doMC)
#' registerDoMC()
#' vim <- varImpact(Y = Y, data = X)
#'
#' ####################################
#' # doSNOW parallel example.
#' library(doSNOW)
#' library(RhpcBLASctl)
#' # Detect the number of physical cores on this computer using RhpcBLASctl.
#' cluster <- makeCluster(get_num_cores())
#' registerDoSNOW(cluster)
#' vim <- varImpact(Y = Y, data = X)
#' stopCluster(cluster)
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
#' doMC::registerDoMC()
#' vim <- varImpact(Y = data$Y, data = subset(data, select=-c(Y, Class, Id)))
#'
#' @export
varImpact = function(Y, data, V = 2,
                     Q.library = c("SL.gam", "SL.glmnet", "SL.mean"),
                     g.library = c("SL.stepAIC"), family = "binomial",
                     minYs = 15, minCell = 0, ncov = 10, corthres = 0.8,
                     impute = "knn",
                     miss.cut = 0.5, verbose=F, parallel = T) {

  # Time the full function execution
  time = system.time({

    # Ensure that Y is numeric; e.g. can't be a factor.
    # (We also assume it's 0/1 but that isn't explictly checked yet.)
    stopifnot(class(Y) %in% c("numeric", "integer"))

    if (family == "binomial" && (min(Y, na.rm = T) < 0 | max(Y, na.rm = T) > 1)) {
      stop("With binomial family Y must be bounded by [0, 1]. Specify family=\"gaussian\" otherwise.")
    }

  ###
  # Get missingness for each column
  # Function for getting total number missing values for vector
  sum.na = function(x) {
    sum(is.na(x))
  }

  ########
  # Applied to Explanatory (X) data frame
  sna = apply(data, 2, sum.na)

  n = nrow(data)

  #######
  # Missing proportion by variable.
  mis.prop = sna / n

  #######
  # Cut-off for eliminating variable for proportion of obs missing.
  data = data[, mis.prop < miss.cut]

  ###### Identify numeric variables (ordered)
  ind.num = sapply(data, is.numeric)

  ## Identify factor variables
  isit.factor = !ind.num

  ###
  # Function that counts # of unique values
  length_unique = function(x) {
    length(unique(x))
  }

  ### num.values is vector of number of unique values by variable
  num.values = apply(data, 2, length_unique)

  ### Function that returns logical T if no variability by variable
  qq.range = function(x, rr) {
    qq = quantile(unclass(x), probs = rr, na.rm = T)
    (qq[2] - qq[1]) == 0
  }

  qq.range.num = function(x, rr) {
    qq = quantile(x, probs = rr, na.rm = T)
    (qq[2] - qq[1]) == 0
  }

  if (sum(isit.factor) == 0) {
    n.fac = 0
  }

  #### data.fac is data frame of variables that are factors
  if (sum(isit.factor) > 0) {
    data.fac = data[, isit.factor, drop = F]

    ## Replace blanks with NA's
    nc = ncol(data.fac)
    for (i in 1:nc) {
      xx = as.character(data.fac[, i])
      xx = factor(xx, exclude = "")
      data.fac[, i] = xx
    }

    # List of column indices to remove.
    dropind = NULL
    for (i in 1:nc) {
      dropind = c(dropind, qq.range(data.fac[, i], rr = c(0.1, 0.9)))
    }

    data.fac = data.fac[, dropind == F, drop = F]

    num.cat = apply(data.fac, 2, length_unique)
    sna = apply(data.fac, 2, sum.na)

    n = nrow(data.fac)
    mis.prop = sna / n

    # POSSIBLE BUG: shouldn't this be using mis.cut rather than 0.5?
    data.fac = data.fac[, mis.prop < 0.5, drop = F]

    # Save how many factors we have in this dataframe.
    num_factors = ncol(data.fac)
  } else {
    num_factors = 0
    datafac.dumW = NULL
  }

  ## Numeric variables
  data.num = data[, !isit.factor, drop = F]

  if (ncol(data.num) > 0) {
    dropind = NULL
    nc = ncol(data.num)
    for (i in 1:nc) {
      # cat(" i = ", i, "\n")
      dropind = c(dropind, qq.range.num(data.num[, i], rr = c(0.1, 0.9)))
    }
    data.num = data.num[, dropind == F, drop = F]
    # Save how many numeric variables we have in this dataframe.
    num_numeric = ncol(data.num)
  } else {
    num_numeric = 0
  }


  if (num_factors > 0) {
    if (verbose) cat("Cleaning up", num_factors, "factor variables.\n")
    ## For each factor, apply qq.range function and get rid of those where
    ## 'true' data.fac is data frame of variables that are factors
    facnames = names(data.fac)
    nam.fac = function(x, name) {
      nc = nchar(x)
      out = paste(name, substr(x, 2, nc), sep = "XX")
      # Remove spaces in variable names.
      out = gsub(" ", "", out)
      return(out)
    }

    newX = NULL
    cumulative_names = NULL

    # TODO: remove this line?
    options(na.action = "na.pass")

    # Loop over each factor variable.
    for (i in 1:num_factors) {
      # cat(' i = ',i,'\n')
      x = data.fac[, i]
      # Convert to a series of indicators.
      # CK: looks like we are omitting the first level?
      inds = model.matrix(~ x - 1)[, -1]
      nmes = colnames(inds)
      if (is.null(nmes)) {
        nmes2 = facnames[i]
      } else {
        # Clean up the names for each indicator.
        nmes2 = nam.fac(nmes, facnames[i])
      }
      # Accumulate the names.
      cumulative_names = c(cumulative_names, nmes2)
      newX = cbind(newX, inds)
    }

    colnames(newX) = cumulative_names

    ## Indexing vector for dummy basis back to original factors
    cc = regexpr("XX", cumulative_names)
    ncc = nchar(cumulative_names)
    cc[cc < 0] = ncc[cc < 0] + 1
    fac.indx = substr(cumulative_names, 1, cc - 1)
    datafac.dum = newX

    ############
    # Missing Basis for Factors
    xp = ncol(datafac.dum)
    sna = apply(datafac.dum, 2, sum.na)
    namesX = colnames(datafac.dum)
    miss.fac = NULL
    names_miss = NULL
    for (k in 1:xp) {
      if (sna[k] > 0) {
        # Again, we are flagging non-missing as 1 and missing as 0 here.
        ix = as.numeric(is.na(datafac.dum[, k]) == F)
        # Replace missing factor indicators with 0.
        datafac.dum[is.na(datafac.dum[, k]), k] = 0
        miss.fac = cbind(miss.fac, ix)
        names_miss = c(names_miss, paste("Imiss_", namesX[k], sep = ""))
      }
    }
    colnames(miss.fac) = names_miss
  } else {
    miss.fac = NULL
    datafac.dum = NULL
  }
  # Finished pre-processing factor variables.
  ###########################################################

  ###########################################################
  # Pre-process numeric/continuous variables.

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
      # TODO: number of bins should be a function argument.
      tst = as.numeric(arules::discretize(Xt, method = "frequency", categories = 10,
                                  ordered = T))
      Xnew = cbind(Xnew, tst)
    }
    colnames(Xnew) = varn
    data.cont.dist = Xnew

    ###############
    # Missing Basis for numeric variables.
    xp = ncol(data.cont.dist)
    n.cont = nrow(data.cont.dist)

    sum.na = function(x) {
      sum(is.na(x))
    }

    sna = apply(data.cont.dist, 2, sum.na)
    nmesX = colnames(data.cont.dist)
    miss.cont = NULL
    nmesm = NULL

    # Create imputed version of the numeric dataframe.
    data.numW = data.num

    # Loop over each binned numeric variable.
    for (k in 1:xp) {
      # Check if that variable has any missing values.
      if (sna[k] > 0) {
        # TODO: is this correct? shouldn't it be is.na() == T?
        # The effect is that the basis is set to 1 if it exists and 0 if it's missing.
        ix = as.numeric(is.na(data.cont.dist[, k]) == F)
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
      impute_info = caret::preProcess(data.num, method=c("medianImpute"))
      data.numW = caret::predict.preProcess(impute_info, data.num)
    } else if (impute == "mean") {
      stop("Mean imputation not implemented yet. Please use another imputation method.")
    } else if (impute == "knn") {
      impute_info = caret::preProcess(data.num, method=c("knnImpute"))
      data.numW = caret::predict.preProcess(impute_info, data.num)
    }

    # Confirm that there are no missing values remaining in data.numW
    stopifnot(sum(is.na(data.numW)) == 0)
  } else {
    impute_info = NULL
    miss.cont = NULL
  }

    cat("Finished pre-processing variables.\n")

    get.tmle.est = function(Y, A, W, delta = NULL, Q.lib, g.lib) {
      ## Because of quirk of program, delete observations with delta=0 if #>0
      ## & < 10
      n = length(Y)
      inc = rep(TRUE, n)
      if (is.null(delta) == F) {
        ss = sum(delta == 0)
        if (ss > 0 & ss < 10) {
          inc[delta == 0] = FALSE
        }
      }
      Y = Y[inc]
      A = A[inc]
      W = W[inc, , drop = F]
      delta = delta[inc]
      tmle.1 = tmle::tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                    Q.SL.library = Q.lib, family = family, verbose = F)
      g1 = tmle.1$g$g1W
      Qst = tmle.1$Qstar[, 2]
      theta = mean(Qst)
      IC = (A/g1) * (Y - Qst) + Qst - theta
      return(list(theta = theta, IC = IC))
    }

    ####
    # Stratified CV to insure balance (by one grouping variable, Y)
    CC.CV = function(V, Y) {
      Ys = unique(Y)
      nys = length(Ys)
      nn = length(Y)
      # Binary outcome so we can do stratified fold generation.
      if (nys == 2) {
        tt = table(Y)
        out = rep(NA, nn)
        for (i in 1:nys) {
          n = as.numeric(tt[i])
          xx = cvTools::cvFolds(n, K = V, R = 1, type = "random")$which
          # TODO: there is a bug in this somewhere.
          out[Y == Ys[i]] = xx
        }
      } else {
        # More than 2 Ys, so don't stratify.
        xx = cvTools::cvFolds(nn, K = V, R = 1, type = "random")$which
        out = xx
      }
      return(out)
    }

    # Create cross-validation folds (2 by default).
    folds = CC.CV(V, Y)

    max.2 = function(x) {
      # Handle missing data manually so we don't get warnings when it occurs.
      x = na.omit(x)
      if (length(x) == 0) {
        -Inf
      } else {
        max(x^2)
      }
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

    n = length(Y)

    ###############################################################################
    # VIM for Factors
    # NOTE: we use && so that conditional will short-circuit if num_factors == 0.
    if (num_factors > 0 && ncol(data.fac) > 0) {
      if (verbose) cat("\nVIM for factors.\n")

      ## Find the level of covariate that has lowest risk
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
      vim_factor = foreach::foreach(i = 1:xc, .verbose=verbose) %do_op% {
      # output <- lapply(1:2, function(i) { pp1=proc.time()
      nameA = names.fac[i]

      if (verbose) cat("i =", i, "Var =", nameA, "out of", xc, "factor variables\n")

      thetaV = NULL
      varICV = NULL
      labV = NULL
      EY0V = NULL
      EY1V = NULL
      nV = NULL

      # Loop over each fold.
      for (kk in 1:V) {
        if (verbose) cat("i =", i, "V =", kk, "\n")

        At = data.fac[folds != kk, i]
        Av = data.fac[folds == kk, i]
        Yt = Y[folds != kk]
        Yv = Y[folds == kk]

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
        W = W[, !apply(is.na(W), 2, all), drop = F]
        Wt = W[folds != kk, , drop = F]
        Wv = W[folds == kk, , drop = F]
        Adum = data.frame(Adum[folds != kk, ])

        ###
        # Pull out any variables that are overly correlated with At (corr coef < corthes)
        #if (sd(Adum) == 0) {
        #  if (verbose) cat("Warning: sd of Adum = 0.\n")
        #}
        corAt = apply(cor(Adum, Wt, use = "complete.obs"), 2, max.2)
        corAt[corAt < -1] = 0
        # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
        incc = abs(corAt) < corthres & is.na(corAt) == F
        Wv = Wv[, incc, drop = F]
        Wt = Wt[, incc, drop = F]

        nw = ncol(Wv)

        ###
        # Skip if number covariates < 10
        if (nw <= 10) {
          Wtsht = Wt
          Wvsht = Wv
        } else {
          if (verbose) cat("Reducing dimensions via clustering.")
          #mydist = as.matrix(hopach::distancematrix(t(Wt), d = "cosangle", na.rm = T))
          mydist = try(hopach::distancematrix(t(Wt), d = "cosangle", na.rm = T),
                       silent = !verbose)
          if (class(mydist) == "try-error") {
            cat("Error in HOPACH clustering: failed to calculate distance matrix.\n")
          }
          hopach.1 = try(hopach::hopach(t(Wt), dmat = mydist, mss = "mean", verbose = F, K = 10,
                                kmax = 3, khigh = 3),
                         silent = !verbose)
          if (class(hopach.1) == "try-error") {
            if (verbose) cat(" Attempt 1 fail.")
            hopach.1 <- try(hopach::hopach(t(Wt), dmat = mydist, mss = "med", verbose = F, K = 10,
                                   kmax = 3, khigh = 3),
                            silent = !verbose)
          }
          if (class(hopach.1) == "try-error") {
            if (verbose) cat(" Attempt 2 fail.")
            #warning("Dimensionality reduction failed. i=", i, "V=", kk, "A=", nameA)
            Wtsht = Wt
            Wvsht = Wv
          } else {
            nlvls = nchar(max(hopach.1$final$labels))
            no = trunc(mean(log10(hopach.1$final$labels)))

            # Find highest level of tree where minimum number of covariates is > ncov
            lvl = 1:nlvls
            ncv = NULL
            for (ii in lvl) {
              ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no - (ii - 1))))))
            }
            ncv = unique(ncv)
            lev = min(min(nlvls, dim(Wt)[2]), min(lvl[ncv >= ncov]))
            two.clust <- unique(trunc(hopach.1$final$labels/(10^(no - (lev - 1)))))
            md <- hopach.1$final$medoids
            mm = md[, 1] %in% two.clust
            incc = md[mm, 2]
            Wtsht = Wt[, incc]
            Wvsht = Wv[, incc]
          }
          if (verbose) cat(" Updated columns:", ncol(Wtsht), "\n")
        }
        # Finished with any needed clustering for variable reduction.

        # cat(' time a = ',proc.time()-pp1,'\n') pp2=proc.time()

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
        minc = apply(table(Av, Yv), 1, min)
        minc2 = apply(table(At, Yt), 1, min)
        if (length(unique(Yt)) == 2) {
          # Binary outcome.
          vals = levA[pmin(minc, minc2) > minCell]
        } else {
          # Continuous outcome.
          vals = levA
        }
        num.cat = length(vals)

        # CK 6/6: don't assume that positive outcome is the rare outcome. (e.g. via table)

        # Number of positive outcomes in training data.
        nYt = sum(Yt[!is.na(At)])
        # Number of positive outcomes in validation data.
        nYv = sum(Yv[!is.na(Av)])

        ############################
        # Don't do if 1) no more than one category of A left or
        # 2) if missingness pattern for A is such that there are few death events left
        # in either (< minYs)
        # Applies only to binary outcomes, not continuous.
        if (length(unique(Yt)) == 2 && (num.cat < 2 || min(nYt, nYv) < minYs)) {
          if (num.cat < 2) {
            error_msg = paste("Skipping", nameA, "due to lack of variation.")
          } else {
            error_msg = paste("Skipping", nameA, "due to minY constraint.", min(nYt, nYv), "<", minYs)
          }
          if (verbose) cat(error_msg, "\n")

          thetaV = c(thetaV, NA)
          varICV = c(varICV, NA)
          labV = rbind(labV, c(NA, NA))
          EY0V = c(EY0V, NA)
          EY1V = c(EY1V, NA)
          nV = c(nV, NA)
        } else {
        # if (num.cat >= 2 & min(nYt, nYv) >= minYs) {
          labmin = NULL
          labmax = NULL
          maxEY1 = -1e+05
          minEY1 = 1e+06
          minj = 0
          maxj = 0
          ICmin = NULL
          ICmax = NULL
          errcnt = 0
          if (verbose) cat("Estimating TMLE on training", paste0("(", num.cat, ")"))
          for (j in 1:num.cat) {
            IA = as.numeric(At == vals[j])
            IA[is.na(IA)] = 0
            # if(min(table(IA,Yt))>=)
            res = try(get.tmle.est(Yt, IA, Wtsht, deltat, Q.lib = Q.library,
                                   g.lib = g.library), silent = TRUE)
            if (class(res) == "try-error") {
              if (verbose) cat("X")
              errcnt = errcnt + 1
            }
            if (class(res) != "try-error") {
              if (verbose) cat(".")
              EY1 = res$theta
              if (EY1 < minEY1) {
                minj = j
                minEY1 = EY1
                labmin = vals[j]
              }
              if (EY1 > maxEY1) {
                maxj = j
                maxEY1 = EY1
                labmax = vals[j]
              }
            }
          }
          if (verbose) cat(" done.\n")
          # cat(' time b = ',proc.time()-pp2,'\n') pp3=proc.time() Now, estimate
          # on validation sample
          if (errcnt == num.cat | minj == maxj) {
            thetaV = c(thetaV, NA)
            varICV = c(varICV, NA)
            labV = rbind(labV, c(NA, NA))
            EY0V = c(EY0V, NA)
            EY1V = c(EY1V, NA)
            nV = c(nV, NA)
          }
          if (errcnt != num.cat & minj != maxj) {
            IA = as.numeric(Av == vals[minj])
            IA[is.na(IA)] = 0
            res = try(get.tmle.est(Yv, IA, Wvsht, deltav, Q.lib = Q.library,
                                   g.lib = c("SL.glmnet", "SL.glm")), silent = TRUE)
            if (class(res) == "try-error") {
              if (verbose) cat("TMLE on validation failed.\n")
              thetaV = c(thetaV, NA)
              varICV = c(varICV, NA)
              labV = rbind(labV, c(NA, NA))
              EY0V = c(EY0V, NA)
              EY1V = c(EY1V, NA)
              nV = c(nV, NA)
            }
            # cat(' time c = ',proc.time()-pp3,'\n') pp4=proc.time()
            if (class(res) != "try-error") {
              IC0 = res$IC
              EY0 = res$theta
              IA = as.numeric(Av == vals[maxj])
              IA[is.na(IA)] = 0
              res2 = try(get.tmle.est(Yv, IA, Wvsht, deltav, Q.lib = Q.library,
                                      g.lib = g.library), silent = TRUE)
              if (class(res2) == "try-error") {
                thetaV = c(thetaV, NA)
                varICV = c(varICV, NA)
                labV = rbind(labV, c(NA, NA))
                EY0V = c(EY0V, NA)
                EY1V = c(EY1V, NA)
                nV = c(nV, NA)
              }
              if (class(res2) != "try-error") {
                IC1 = res2$IC
                EY1 = res2$theta
                thetaV = c(thetaV, EY1 - EY0)
                varICV = c(varICV, var(IC1 - IC0))
                labV = rbind(labV, c(labmin, labmax))
                EY0V = c(EY0V, EY0)
                EY1V = c(EY1V, EY1)
                nV = c(nV, length(Yv))
              }
            }
          }
        }

      }

      list(EY1V, EY0V, thetaV, varICV, labV, nV, "factor")
      # print(data.frame(EY1V,EY0V,thetaV,varICV,labV,nV))
    } # End foreach loop.

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
      (cor(na.omit(cbind(x, y)))[1, 2])^2
    }

    ######################################################
    # We use && so that the second check will be skipped when num_numeric == 0.
    if (num_numeric > 0 && ncol(data.cont.dist) > 0) {
      if (verbose) cat("\nVIM for numeric variables.\n")

      xc = ncol(data.cont.dist)
      names.cont = colnames(data.cont.dist)
      n.cont = nrow(data.cont.dist)

      numcat.cont = apply(data.cont.dist, 2, length_unique)
      xc = length(numcat.cont)

      cats.cont = lapply(1:xc, function(i) {
        sort(unique(data.cont.dist[, i]))
      })

      ### Loop over each numeric variable.
      #vim_numeric = lapply(1:xc, function(i) {
      vim_numeric = foreach::foreach(i = 1:xc) %do_op% {
      nameA = names.cont[i]

      if (verbose) cat("i =", i, "Var =", nameA, "out of", xc, "numeric variables\n")

      thetaV = NULL
      varICV = NULL
      labV = NULL
      EY0V = NULL
      EY1V = NULL
      nV = NULL

      for (kk in 1:V) {
        if (verbose) cat("v =", kk, "out of V =", V, "\n")

        At = data.cont.dist[folds != kk, i]
        Av = data.cont.dist[folds == kk, i]
        Yt = Y[folds != kk]
        Yv = Y[folds == kk]

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
        if (singleAY1 || length(na.omit(unique(Atnew))) <= 1 || length(na.omit(unique(Avnew))) <= 1) {
          error_msg = paste("Skipping", nameA, "in this fold because there is no variation.")
          if (verbose) cat(error_msg, "\n")
          warning(error_msg)
          thetaV = c(thetaV, NA)
          varICV = c(varICV, NA)
          labV = rbind(labV, c(NA, NA))
          EY0V = c(EY0V, NA)
          EY1V = c(EY1V, NA)
          nV = c(nV, NA)
        } else {
        #if (length(na.omit(unique(Atnew))) > 1 & length(na.omit(unique(Avnew))) > 1) {
          labs = names(table(Atnew))
          Atnew = as.numeric(Atnew) - 1
          Avnew = as.numeric(Avnew) - 1
          numcat.cont[i] = length(labs)
          # change this to match what was done for factors - once
          # cats.cont[[i]]=as.numeric(na.omit(unique(Atnew)))
          cats.cont[[i]] = as.numeric(names(table(Atnew)))
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
          W = data.frame(data.numW[, -i, drop = F], miss.cont, datafac.dumW, miss.fac)
          W = W[, !apply(is.na(W), 2, all), drop = F]
          Wt = W[folds != kk, , drop = F]
          Wv = W[folds == kk, , drop = F]

          nmesW = names(Wt)
          mtch = match(nmesW, paste0("Imiss_", nameA))
          Wt = Wt[, is.na(mtch), drop = F]
          Wv = Wv[, is.na(mtch), drop = F]
          ### Pull out any variables that are overly correlated with At (corr coef
          ### < corthes)
          corAt = apply(Wt, 2, cor.two, y = At)
          # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
          incc = corAt < corthres & is.na(corAt) == F
          Wv = Wv[, incc, drop = F]
          Wt = Wt[, incc, drop = F]
          #### Use HOPACH to reduce dimension of W to some level of tree
          nw = ncol(Wv)

          ### Skip if number covariates < 10
          if (nw <= 10) {
            Wtsht = Wt
            Wvsht = Wv
          } else {
            if (verbose) cat("Reducing dimensions via clustering.")
          #if (nw > 10) {
            #mydist = as.matrix(hopach::distancematrix(t(Wt), d = "cosangle", na.rm = T))
            mydist = try(hopach::distancematrix(t(Wt), d = "cosangle", na.rm = T),
                         silent = !verbose)
            if (class(mydist) == "try-error") {
              cat("Error in HOPACH clustering: failed to calculate distance matrix.\n")
            }
            hopach.1 = try(hopach::hopach(t(Wt), dmat = mydist, mss = "mean", verbose = F),
                           silent = !verbose)
            if (class(hopach.1) == "try-error") {
              if (verbose) cat(" Attempt 1 fail.")
              # Retry with a different specification.
              hopach.1 <- try(hopach::hopach(t(Wt), dmat = mydist, mss = "med", verbose = F),
                              silent = !verbose)
            }
            if (class(hopach.1) == "try-error") {
              if (verbose) cat(" Attempt 2 fail.")
              Wtsht = Wt
              Wvsht = Wv
            } else {
            #if (class(hopach.1) != "try-error") {
              nlvls = nchar(max(hopach.1$final$labels))
              no = trunc(mean(log10(hopach.1$final$labels)))
              # Find highest level of tree where minimum number of covariates is > ncov
              lvl = 1:nlvls
              ncv = NULL
              for (ii in lvl) {
                ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no - (ii - 1))))))
              }
              ncv = unique(ncv)
              lev = min(min(nlvls, dim(Wt)[2]), min(lvl[ncv >= ncov]))
              two.clust <- unique(trunc(hopach.1$final$labels/(10^(no - (lev - 1)))))
              md = hopach.1$final$medoids
              mm = md[, 1] %in% two.clust
              incc = md[mm, 2]
              Wtsht = Wt[, incc]
              Wvsht = Wv[, incc]
            }
            if (verbose) cat(" Updated columns:", ncol(Wtsht), "\n")
          }
          deltat = as.numeric(is.na(Yt) == F & is.na(Atnew) == F)
          deltav = as.numeric(is.na(Yv) == F & is.na(Avnew) == F)

          if (sum(deltat == 0) < 10) {
            Yt = Yt[deltat == 1]
            Wtsht = Wtsht[deltat == 1, ]
            Atnew = Atnew[deltat == 1]
            deltat = deltat[deltat == 1]
          }

          vals = cats.cont[[i]]
          maxEY1 = -1e+05
          minEY1 = 1e+06
          minj = 0
          maxj = 0
          ICmin = NULL
          ICmax = NULL
          Atnew[is.na(Atnew)] = -1
          Avnew[is.na(Avnew)] = -1
          # Only applies to binary outcomes.
          if (length(unique(Yt)) == 2 && min(table(Avnew[Avnew >= 0], Yv[Avnew >= 0])) <= minCell) {
            error_msg = paste("Skipping", nameA, "due to minCell constraint.\n")
            if (verbose) cat(error_msg)
            warning(error_msg)
            thetaV = c(thetaV, NA)
            varICV = c(varICV, NA)
            labV = rbind(labV, c(NA, NA))
            EY0V = c(EY0V, NA)
            EY1V = c(EY1V, NA)
            nV = c(nV, NA)
          }
          # CK TODO: this is not exactly the opposite of the IF above. Is that intentional?
          if (length(unique(Yt)) > 2 || min(table(Avnew, Yv)) > minCell) {
            labmin = NULL
            labmax = NULL
            errcnt = 0
            if (verbose) cat("Estimating training TMLEs", paste0("(", numcat.cont[i], ")"))
            for (j in 1:numcat.cont[i]) {
              # cat(' i = ',i,' kk = ',kk,' j = ',j,'\n')
              IA = as.numeric(Atnew == vals[j])
              res = try(get.tmle.est(Yt, IA, Wtsht, deltat, Q.lib = Q.library,
                                     g.lib = g.library), silent = T)
              if (class(res) == "try-error") {
                # Error.
                if (verbose) cat("X")
                errcnt = errcnt + 1
              } else {
                if (verbose) cat(".")
              #if (class(res) != "try-error") {
                EY1 = res$theta
                if (EY1 < minEY1) {
                  minj = j
                  minEY1 = EY1
                  labmin = labs[j]
                }
                # CK 6/6: convert to else? But what about equality?
                if (EY1 > maxEY1) {
                  maxj = j
                  maxEY1 = EY1
                  labmax = labs[j]
                }
              }
            }
            if (verbose) cat("done.\n")

            #####
            # Now, estimate on validation sample
            if (errcnt == numcat.cont[i] | minj == maxj) {
              if (verbose) {
                if (minj == maxj) {
                  cat("Skipping", nameA, "because min level = max level =", minj, "\n")
                } else {
                  cat("Skipping", nameA, "because no level succeeded\n")
                }
              }
              thetaV = c(thetaV, NA)
              varICV = c(varICV, NA)
              labV = rbind(labV, c(NA, NA))
              EY0V = c(EY0V, NA)
              EY1V = c(EY1V, NA)
              nV = c(nV, NA)
            }
            if (errcnt != numcat.cont[i] & minj != maxj) {
              # Estimate with the minimum level.
              IA = as.numeric(Avnew == vals[minj])
              if (verbose) cat("Estimate on validation: ")
              res = try(get.tmle.est(Yv, IA, Wvsht, deltav, Q.lib = Q.library,
                                     g.lib = g.library), silent = T)
              if (verbose) cat("min ")
              if (class(res) == "try-error") {
                if (verbose) cat(" Failed :/\n")
                thetaV = c(thetaV, NA)
                varICV = c(varICV, NA)
                labV = rbind(labV, c(NA, NA))
                EY0V = c(EY0V, NA)
                EY1V = c(EY1V, NA)
                nV = c(nV, NA)
              }
              if (class(res) != "try-error") {
                # Estimate with the maximum level.
                IC0 = res$IC
                EY0 = res$theta
                IA = as.numeric(Avnew == vals[maxj])
                res2 = try(get.tmle.est(Yv, IA, Wvsht, deltav,
                                        Q.lib = Q.library, g.lib = g.library), silent = TRUE)
                if (verbose) cat("max")
                if (class(res2) == "try-error") {
                  if (verbose) cat(" Failed :/\n")
                  thetaV = c(thetaV, NA)
                  varICV = c(varICV, NA)
                  labV = rbind(labV, c(NA, NA))
                  EY0V = c(EY0V, NA)
                  EY1V = c(EY1V, NA)
                  nV = c(nV, NA)
                }
                if (class(res2) != "try-error") {
                  if (verbose) cat(".\n")
                  IC1 = res2$IC
                  EY1 = res2$theta
                  thetaV = c(thetaV, EY1 - EY0)
                  varICV = c(varICV, var(IC1 - IC0))
                  labV = rbind(labV, c(labmin, labmax))
                  EY0V = c(EY0V, EY0)
                  EY1V = c(EY1V, EY1)
                  nV = c(nV, length(Yv))
                }
              }
            }
          }
        }
      }
      list(EY1V, EY0V, thetaV, varICV, labV, nV, "numeric")
    } # end foreach loop.

      if (verbose) cat("Numeric VIMs:", length(vim_numeric), "\n")

      # Confirm that we have the correct number of results, otherwise fail out.
      stopifnot(length(vim_numeric) == xc)

      colnames_numeric = colnames(data.cont.dist)
    } else {
      colnames_numeric = NULL
      vim_numeric = NULL
      data.numW = NULL
      cat("No numeric variables - skim VIM estimation.\n")
    }

    if (verbose) cat("Completed VIM estimation.\n")

    # Combine the separate continuous and factor results.

    num_numeric = length(colnames_numeric)
    num_factor = length(colnames_factor)

    factr = c(rep("ordered", num_numeric), rep("factor", num_factor))
    vars = c(colnames_numeric, colnames_factor)

    vim_combined = c(vim_numeric, vim_factor)
    names(vim_combined) = vars

    lngth = sapply(vim_combined, function(x) length(x))

    # Restrict to vim results that have 7 elements.
    vim_combined = vim_combined[lngth == 7]
    vars = vars[lngth == 7]
    factr = factr[lngth == 7]

    # Set defaults for variables we want to return.
    outres.all = outres.cons = outres.byV = NULL

    # Get rid of any variables that have a validation sample with no
    # estimates of variable importance.
    if (length(vim_combined) == 0) {
      error_msg = "No VIMs could be calculated due to sample size, etc."
      if (verbose) cat(error_msg, "\n")
      warning(error_msg)
      # TODO: also write output to the file in a separate function call.
      #write("No VIM's could be calculated due to sample size, etc",
      #     file = "AllReslts.tex")
    } else {
      lngth2 = sapply(vim_combined, function(x) length(na.omit(x[[1]])))
      out.put = vim_combined[lngth2 == V]
      if (length(out.put) == 0) {
        error_msg = "No VIMs could be calculated due to sample size, etc."
        if (verbose) cat(error_msg, "\n")
        warning(error_msg)
      } else {
        vars = vars[lngth2 == V]
        factr = factr[lngth2 == V]
        tst = lapply(out.put, function(x) x[[3]])
        tst = do.call(rbind, tst)
        tot.na = function(x) {
          sum(is.na(x))
        }
        if (verbose) cat("Dim tst:", paste(dim(tst)), "\n")
        xx = apply(tst, 1, tot.na)
        out.sht = out.put[xx == 0]
        vars = vars[xx == 0]
        factr = factr[xx == 0]

        # names(out.sht)=vars[xx==0]
        tst = lapply(out.sht, function(x) x[[1]])
        EY1 = do.call(rbind, tst)
        tst = lapply(out.sht, function(x) x[[2]])
        EY0 = do.call(rbind, tst)
        tst = lapply(out.sht, function(x) x[[3]])
        theta = do.call(rbind, tst)
        tst = lapply(out.sht, function(x) x[[4]])
        varIC = do.call(rbind, tst)
        tst = lapply(out.sht, function(x) x[[6]])
        nV = do.call(rbind, tst)
        n = sum(nV[1, ])
        SEV = sqrt(varIC/nV)
        ##### Get labels for each of the training sample
        labs.get = function(x, fold) {
          lbel = rep(1:fold, 2)
          oo = order(lbel)
          lbel = lbel[oo]
          out = as.vector(t(x))
          names(out) = paste("v.", lbel, rep(c("a_L", "a_H"), 2), sep = "")
          out
        }
        tst = lapply(out.sht, function(x) x[[5]])
        tst = lapply(tst, labs.get, fold = V)
        lbs = do.call(rbind, tst)

        meanvarIC = apply(varIC, 1, mean)
        psi = apply(theta, 1, mean)
        SE = sqrt(meanvarIC/n)

        lower = psi - 1.96 * SE
        upper = psi + 1.96 * SE
        signdig = 3
        CI95 = paste0("(", signif(lower, signdig), " - ", signif(upper, signdig), ")")

        # 1-sided p-value
        pvalue = 1 - pnorm(psi/SE)
        ##### FOR THETA (generalize to chi-square test?)  TT =
        ##### (theta[,1]-theta[,2])/sqrt(SEV[,1]^2+SEV[,2]^2)
        ##### pval.comp=2*(1-pnorm(abs(TT))) FOR levels (just make sure in same
        ##### order)
        nc = sum(factr == "ordered")
        n = length(factr)
        ### Ordered variables first
        length.uniq = function(x) {
          length(unique(x))
        }
        cons = NULL
        if (nc > 0) {
          dir = NULL
          for (i in 1:V) {
            ltemp = lbs[1:nc, i * 2 - 1]
            xx = regexpr(",", ltemp)
            lwr = as.numeric(substr(ltemp, 2, xx - 1))
            utemp = lbs[1:nc, i * 2]
            xx = regexpr(",", utemp)
            nx = nchar(utemp)
            uwr = as.numeric(substr(utemp, xx + 1, nx - 1))
            dir = cbind(dir, uwr > lwr)
          }
          cons = apply(dir, 1, length.uniq)
        }
        #### Factors
        if (n - nc > 0) {
          lwr = NULL
          uwr = NULL
          for (i in 1:V) {
            lwr = cbind(lwr, lbs[(nc + 1):n, i * 2 - 1])
            uwr = cbind(uwr, lbs[(nc + 1):n, i * 2 - 1])
          }
          conslwr = apply(lwr, 1, length.uniq)
          consupr = apply(uwr, 1, length.uniq)
          cons = c(cons, conslwr * consupr)
        }
        # consist= (cons==1 & pval.comp > 0.05)
        signsum = function(x) {
          sum(sign(x))
        }
        consist = cons == 1 & abs(apply(theta, 1, signsum)) == V
        procs = c("Holm", "BH")
        if (n > 1) {
          res = multtest::mt.rawp2adjp(pvalue, procs)
          oo = res$index
          outres = data.frame(factor = factr[oo], theta[oo, ],
                              psi[oo], CI95[oo], res$adj, lbs[oo, ], consist[oo])
        }
        if (n == 1) {
          outres = data.frame(factor = factr, theta, psi, CI95,
                              rawp = pvalue, Holm = pvalue, BH = pvalue, lbs, consist)
        }

        # Restrict to variables that aren't missing their p-value.
        outres = outres[!is.na(outres[, "rawp"]), , drop = F]

        names(outres)[1:(1 + 2 * V)] = c("VarType", paste0("psiV", 1:V), "AvePsi", "CI95")
        names(outres)[(9 + 2 * V)] = "Consistent"

        ### Get Consistency Measure and only significant Make BH cut-off flexible
        ### in future versions (default at 0.05)
        outres.cons = outres[outres[, "BH"] < 0.05, , drop = F]
        outres.cons = outres.cons[outres.cons[, "Consistent"],
                                  c("VarType", "AvePsi", "rawp"), drop = F]
        colnames(outres.cons) = c("Type", "Estimate", "P-value")


        # drops = c('VarType','description','Holm,')
        # outres.all=outres[,!(names(outres) %in% drops)]
        outres.byV = outres[, c(2:(2 + V - 1), 9:(9 + 2 * V)), drop = F]
        outres.all = outres[, c("VarType", "AvePsi", "CI95", "rawp", "BH"), drop = F]
        colnames(outres.all) = c("Type", "Estimate", "CI95", "P-value", "Adj. p-value")
      }
    }

  }) # End timing the full execution.

  # Return results.
  results = list(results_consistent = outres.cons,
                 results_all = outres.all,
                 results_by_fold = outres.byV,
                 V=V, g.library = g.library, Q.library = Q.library,
                 minCell = minCell, minYs = minYs,
                 family=family,  datafac.dumW  = datafac.dumW,
                 miss.fac = miss.fac,
                 data.numW = data.numW, impute_info = impute_info,
                 time = time)
  # Set a custom class so that we can override print and summary.
  class(results) = "varImpact"
  invisible(results)
}
