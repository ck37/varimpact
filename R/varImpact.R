#' Variable Impact Estimation
#' \code{varImpact} returns variable importance statistics ordered
#' by statistical significance
#'
#' Returns ordered estimates of variable importance measures using
#' a combination of data-adaptive target parameter estimation, and
#' targeted maximum likelihood estimation (TMLE).
#'
#'  The function performs the following functions.
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
#'  \item Returns 3 latex table files:
#'  \itemize{
#'    \item AllReslts.tex - the file with cross-validated average variable impacts ordered by statistical significance.
#'    \item byV.tex - the comparison levels used within each validation sample.  Either integer ordering of factors or short-hand for percentile cut-off (0-1 is the 10th percentile, 10+ is the 100th percentile)
#'    \item ConsistReslts.tex - the ``consistent'' significant results, meaning those with consistent categories chosen as comparison groups among factors and consistent ordering for numeric variables.
#' }
#'  \item Things to do include making options to avoid errors include putting
#' minimum cell size on validation sample of A vs. Y
#' and implementing CV-TMLE (minCell), making examples, putting authors
#' and references and see also's.  Allow reporting of results that
#' randomly do not have estimates for some of validation samples.
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
#' @param cothres cut-off correlation with explanatory
#' variable for inclusion of an adjustment variables
#' @param   dirout directory to write output
#' @param   outname prefix name for files
#' @param miss.cut eliminates explanatory (X) variables with proportion
#' of missing obs > cut.off

varImpact = function(Y, data, V,
                     Q.library = c("SL.gam", "SL.glmnet", "SL.mean"),
                     g.library = c("SL.stepAIC"), family = "binomial",
                     minYs = 15, minCell = 0, ncov = 10, corthres = 0.8,
                     dirout = NULL, outname = NULL,
                     miss.cut = 0.5) {

  # TODO: remove this line and replace all "fam" referneces with "family".
  fam = family
  # TODO: remove this line and replace all "data1" references with "data".
  data1 = data

  ###################################
  # OUTSTANDING PROGRAMMING ISSUES (CK: possibly old)
  # 1) How to add 'null' columns to data frame?

  ###
  # Get missingness for each column
  # Function for getting total number missing values for vector
  sum.na = function(x) {
    sum(is.na(x))
  }

  ########
  # Applied to Explanatory (X) data frame
  sna = apply(data1, 2, sum.na)

  ####### n = # of row
  n = dim(data1)[1]

  #######
  # Missing proportion by variable.
  mis.prop = sna / n

  #######
  # Cut-off for eliminating variable for proportion of obs missing.
  data1 = data1[, mis.prop < miss.cut]

  ###### Identify numeric variables (ordered)
  ind.num = sapply(data1, is.numeric)

  ## Identify factor variables
  # isit.factor = sapply(data1,is.factor)
  isit.factor = !ind.num

  # Number of columns in the reduced dataframe.
  ndata1 = dim(data1)[2]

  ###
  # Function that counts # of unique values
  num.val = function(x) {
    length(unique(x))
  }
  ### num.values is vector of number of unique values by variable
  num.values = apply(data1, 2, num.val)

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
    data.fac = data1[, isit.factor, drop = F]

    ## Replace blanks with NA's
    nc = dim(data.fac)[2]
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

    # TODO: remove this function, as it duplicates num.val() defined above.
    ln.unique = function(x) {
      length(unique(x))
    }
    num.cat = apply(data.fac, 2, ln.unique)
    sna = apply(data.fac, 2, sum.na)

    n = dim(data.fac)[1]
    mis.prop = sna / n

    # POSSIBLE BUG: shouldn't this be using mis.cut rather than 0.5?
    data.fac = data.fac[, mis.prop < 0.5, drop = F]

    # Save how many factors we have in this dataframe.
    n.fac = dim(data.fac)[2]
  }

  ## Numeric variables
  data.num = data1[, isit.factor == F, drop = F]
  if ((dim(data.num)[2]) == 0) {
    n.num = 0
  }

  if ((dim(data.num)[2]) > 0) {
    dropind = NULL
    nc = dim(data.num)[2]
    for (i in 1:nc) {
      cat(" i = ", i, "\n")
      dropind = c(dropind, qq.range.num(data.num[, i], rr = c(0.1, 0.9)))
    }
    data.num = data.num[, dropind == F, drop = F]
    n.num = dim(data.num)[2]
  }

  ###########################################################################
  # Do following duplicate code modified if there are no factors
  # Need to re-code to remove redundancy
  ################################################

  if (n.num > 0 & n.fac == 0) {
    # TODO: remove this, already have num.val() above.
    ln.unique = function(x) {
      length(unique(x))
    }
    X = data.num
    xc = dim(X)[2]
    qt = apply(na.omit(X), 2, quantile, probs = seq(0.1, 0.9, 0.1))

    newX = NULL
    coln = NULL

    varn = colnames(X)

    num.cat = apply(X, 2, ln.unique)

    ### Processing continuous variables
    Xnew = NULL
    for (k in 1:xc) {
      Xt = X[, k]
      tst = as.numeric(arules::discretize(Xt, method = "frequency", categories = 10, ordered = T))
      Xnew = cbind(Xnew, tst)
    }
    colnames(Xnew) = varn
    data.cont.dist = Xnew

    ###############
    # Missing Basis
    xp = dim(data.cont.dist)[2]
    sum.na = function(x) {
      sum(is.na(x))
    }
    sna = apply(data.cont.dist, 2, sum.na)
    nmesX = colnames(data.cont.dist)

    miss.cont = NULL
    nmesm = NULL

    for (k in 1:xp) {
      if (sna[k] > 0) {
        ix = as.numeric(is.na(data.cont.dist[, k]) == F)
        miss.cont = cbind(miss.cont, ix)
        nmesm = c(nmesm, paste("Imiss_", nmesX[k], sep = ""))
      }
    }
    colnames(miss.cont) = nmesm

    numcat.cont = apply(data.cont.dist, 2, ln.unique)

    xc = length(numcat.cont)
    cats.cont <- lapply(1:xc, function(i) {
      sort(unique(data.cont.dist[, i]))
    })


    ## Find the level of covariate that has lowest risk
    data.numW = data.num
    data.numW[is.na(data.num)] = 0
    n = length(Y)
    get.tmle.est = function(Y, A, W, delta = NULL, Q.lib, g.lib) {
      # Because of quirk of program, delete obs with delta=0 if #>0  & < 10
      n = length(Y)

      # Vector of whether or not to include an obseration.
      inc = rep(TRUE, n)
      if (is.null(delta) == F) {
        ss = sum(delta == 0)
        if (ss > 0 & ss < 10) {
          inc[delta == 0] = FALSE
        }
      }

      # Subset to our desired rows.
      Y = Y[inc]
      A = A[inc]
      W = W[inc, , drop = F]
      delta = delta[inc]

      tmle.1 = tmle::tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                    Q.SL.library = Q.lib, family = fam, verbose = FALSE)
      g1 = tmle.1$g$g1W
      Qst = tmle.1$Qstar[, 2]
      theta = mean(Qst)
      IC = (A/g1) * (Y - Qst) + Qst - theta
      # TODO: return other aspects of the TMLE results, e.g. g1 or Qstar?
      return(list(theta = theta, IC = IC))
    }

    #### Stratified CV to insure balance (by one grouping variable, Y)
    CC.CV = function(V, Y) {
      nn = length(Y)
      tt = table(Y)
      Ys = unique(Y)
      nys = length(Ys)
      out = rep(NA, nn)
      for (i in 1:nys) {
        n = as.numeric(tt[i])
        xx = cvTools::cvFolds(n, K = V, R = 1, type = "random")$which
        out[Y == Ys[i]] = xx
      }
      return(out)
    }

    folds = CC.CV(V, Y)

    max.2 = function(x) {
      max(x^2, na.rm = T)
    }

    ###############################################################################
    cor.two = function(x, y) {
      (cor(na.omit(cbind(x, y)))[1, 2])^2
    }

    # TODO: delete this line? We already made the folds 10 lines up.
    folds = CC.CV(V, Y)

    # detach('package:cvTools', unload=TRUE)
    # detach('package:lattice', unload=TRUE)
    # library(doParallel)
    # cl <- makeCluster(20)
    # registerDoParallel(cl)

    names.cont = colnames(data.cont.dist)
    xc = dim(data.cont.dist)[2]
    n.cont = dim(data.cont.dist)[1]

    ###########################################################################

    # Loop over each column in X to conduct the VIM analysis.
    out.put <- lapply(1:xc, function(i) {
      nameA = names.cont[i]
      thetaV = NULL
      varICV = NULL
      labV = NULL
      EY0V = NULL
      EY1V = NULL
      nV = NULL

      # Loop over cross-validation folds.
      for (kk in 1:V) {
        # cat(' i = ',i,' V = ',kk,'\n')

        # A in the training set.
        At = data.cont.dist[folds != kk, i]

        # A in the validation set.
        Av = data.cont.dist[folds == kk, i]

        # Y in the training set.
        Yt = Y[folds != kk]

        # Y in the validation set.
        Yv = Y[folds == kk]

        AY1 = At[Yt == 1 & is.na(At) == F]

        # Create a histogram with automatic bin selection.
        hh = histogram::histogram(AY1, verbose = F, type = "irregular")$breaks

        # Ensure that the min and max breaks bound the distribution of At.
        if (hh[length(hh)] < max(At, na.rm = T)) {
          hh[length(hh)] = max(At, na.rm = T) + 0.1
        }
        if (hh[1] > min(At[At > 0], na.rm = T)) {
          hh[1] = min(At[At > 0], na.rm = T) - 0.1
        }

        # Discretize A using the histogram breaks.
        Atnew = cut(At, breaks = hh)
        Avnew = cut(Av, breaks = hh)

        if (length(na.omit(unique(Atnew))) <= 1 | length(na.omit(unique(Avnew))) <= 1) {
          thetaV = c(thetaV, NA)
          varICV = c(varICV, NA)
          labV = rbind(labV, c(NA, NA))
          EY0V = c(EY0V, NA)
          EY1V = c(EY1V, NA)
          nV = c(nV, NA)
        }
        if (length(na.omit(unique(Atnew))) > 1 & length(na.omit(unique(Avnew))) > 1) {
          labs = names(table(Atnew))
          Atnew = as.numeric(Atnew) - 1
          Avnew = as.numeric(Avnew) - 1
          numcat.cont[i] = length(labs)
          # change this to match what was done for factors - once
          # cats.cont[[i]]=as.numeric(na.omit(unique(Atnew)))
          cats.cont[[i]] = as.numeric(names(table(Atnew)))
          ### acit.numW is just same as data.cont.dist except with NA's replaced by
          ### 0's.
          if (is.null(miss.cont)) {
            miss.cont = matrix(0, n.cont, 1)
          }

          # W restricted to the training set.
          Wt = data.frame(data.numW[folds != kk, -i], miss.cont[folds != kk, ])
          nmesW = names(Wt)

          # TODO: convert to paste0
          mtch = match(nmesW, paste("Imiss_", nameA, sep = ""))

          Wt = Wt[, is.na(mtch)]
          Wv = data.frame(data.numW[folds == kk, -i], miss.cont[folds == kk, ])
          Wv = Wv[, is.na(mtch)]

          ###
          # Pull out any variables that are overly correlated with At (corr coef < corthes)
          corAt = apply(Wt, 2, cor.two, y = At)
          # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
          incc = corAt < corthres & is.na(corAt) == F

          Wv = Wv[, incc]
          Wt = Wt[, incc]

          # ERROR: what is nw referring to? data.numW?
          # CK 5/31: let's guess at what nw should be.
          nw = ncol(Wt)

          # TODO: this should be a user-configurable parameter.
          if (nw <= 10) {
            Wtsht = Wt
            Wvsht = Wv
          }
          if (nw > 10) {
            mydist <- as.matrix(distancematrix(t(Wt), d = "cosangle", na.rm = T))
            hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "mean",
                                   verbose = FALSE), silent = TRUE)
            if (class(hopach.1) == "try-error") {
              hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "med",
                                     verbose = FALSE), silent = TRUE)
            }
            if (class(hopach.1) == "try-error") {
              Wtsht = Wt
              Wvsht = Wv
            }
            if (class(hopach.1) != "try-error") {
              nlvls = nchar(max(hopach.1$final$labels))
              no <- trunc(mean(log10(hopach.1$final$labels)))
              ### Find highest level of tree where minimum number of covariates is >
              ### ncov
              lvl = 1:nlvls
              ncv = NULL
              for (ii in lvl) {
                ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no -
                                                                             (ii - 1))))))
              }
              ncv = unique(ncv)
              lev = min(min(nlvls, dim(Wt)[2]), min(lvl[ncv >=
                                                          ncov]))
              two.clust <- unique(trunc(hopach.1$final$labels/(10^(no -
                                                                     (lev - 1)))))
              md <- hopach.1$final$medoids
              mm = md[, 1] %in% two.clust
              incc = md[mm, 2]
              Wtsht = Wt[, incc]
              Wvsht = Wv[, incc]
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

          vals = cats.cont[[i]]
          maxEY1 = -1e+05
          minEY1 = 1e+06
          minj = 0
          maxj = 0
          ICmin = NULL
          ICmax = NULL
          Atnew[is.na(Atnew)] = -1
          Avnew[is.na(Avnew)] = -1
          if (min(table(Avnew[Avnew >= 0], Yv[Avnew >= 0])) <= minCell) {
            thetaV = c(thetaV, NA)
            varICV = c(varICV, NA)
            labV = rbind(labV, c(NA, NA))
            EY0V = c(EY0V, NA)
            EY1V = c(EY1V, NA)
            nV = c(nV, NA)
          }
          if (min(table(Avnew, Yv)) > minCell) {
            labmin = NULL
            labmax = NULL
            errcnt = 0
            for (j in 1:numcat.cont[i]) {
              # cat(' i = ',i,' kk = ',kk,' j = ',j,'\n')
              IA = as.numeric(Atnew == vals[j])
              res = try(get.tmle.est(Yt, IA, Wtsht, deltat, Q.lib = Q.library,
                                     g.lib = g.library), silent = T)
              if (class(res) == "try-error") {
                errcnt = errcnt + 1
              }
              if (class(res) != "try-error") {
                EY1 = res$theta
                if (EY1 < minEY1) {
                  minj = j
                  minEY1 = EY1
                  labmin = labs[j]
                }
                if (EY1 > maxEY1) {
                  maxj = j
                  maxEY1 = EY1
                  labmax = labs[j]
                }
              }
            }
            ##### Now, estimate on validation sample
            if (errcnt == numcat.cont[i] | minj == maxj) {
              thetaV = c(thetaV, NA)
              varICV = c(varICV, NA)
              labV = rbind(labV, c(NA, NA))
              EY0V = c(EY0V, NA)
              EY1V = c(EY1V, NA)
              nV = c(nV, NA)
            }
            if (errcnt != numcat.cont[i] & minj != maxj) {
              IA = as.numeric(Avnew == vals[minj])
              res = try(get.tmle.est(Yv, IA, Wvsht, deltav, Q.lib = Q.library,
                                     g.lib = g.library), silent = T)
              if (class(res) == "try-error") {
                thetaV = c(thetaV, NA)
                varICV = c(varICV, NA)
                labV = rbind(labV, c(NA, NA))
                EY0V = c(EY0V, NA)
                EY1V = c(EY1V, NA)
                nV = c(nV, NA)
              }
              if (class(res) != "try-error") {
                IC0 = res$IC
                EY0 = res$theta
                IA = as.numeric(Avnew == vals[maxj])
                res2 = try(get.tmle.est(Yv, IA, Wvsht, deltav,
                                        Q.lib = Q.library, g.lib = g.library), silent = T)
                if (class(res2) == "try-error") {
                  thetaV = c(thetaV, NA)
                  varICV = c(varICV, NA)
                  labV = rbind(labV, c(NA, NA))
                  EY0V = c(EY0V, NA)
                  EY1V = c(EY1V, NA)
                  nV = c(nV, NA)
                }
                if (class(res2) != "try-error") {
                  IC1 = res$IC
                  EY1 = res$theta
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
      # print(list(EY1V,EY0V,thetaV,varICV,labV,nV))
      list(EY1V, EY0V, thetaV, varICV, labV, nV)
    })

    vars.cont = colnames(data.cont.dist)
    stf = c("vars.cont", "out.cont")

    # vars=colnames(data.fac)
    nc = length(vars.cont)
    # nf=length(vars)
    factr = rep("ordered", nc)
    vars = vars.cont
    names(out.put) = vars
    lngth = sapply(out.put, function(x) length(x))
    out.put = out.put[lngth == 6]
    vars = vars[lngth == 6]
    factr = factr[lngth == 6]
    lngth2 = sapply(out.put, function(x) length(na.omit(x[[1]])))
    out.put = out.put[lngth2 == V]
    vars = vars[lngth2 == V]
    factr = factr[lngth2 == V]

    tst = lapply(out.put, function(x) x[[3]])
    tst = do.call(rbind, tst)
    tot.na = function(x) {
      sum(is.na(x))
    }
    xx = apply(tst, 1, tot.na)
    out.sht = out.put[xx == 0]
    vars = vars[xx == 0]
    factr = factr[xx == 0]

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
    tst = lapply(tst, labs.get, fold = 2)
    lbs = do.call(rbind, tst)


    meanvarIC = apply(varIC, 1, mean)
    psi = apply(theta, 1, mean)
    SE = sqrt(meanvarIC/n)
    lower = psi - 1.96 * SE
    upper = psi + 1.96 * SE
    signdig = 3
    CI95 = paste("(", signif(lower, signdig), " - ", signif(upper,
                                                            signdig), ")", sep = "")
    # 1-sided p-value
    pvalue = 1 - pnorm(psi/SE)

    #####
    # FOR THETA (generalize to chi-square test?)
    # TT = (theta[,1]-theta[,2])/sqrt(SEV[,1]^2+SEV[,2]^2)
    # pval.comp = 2*(1-pnorm(abs(TT)))


    #####
    # FOR levels (just make sure in same order)
    nc = n = length(factr)
    dir = NULL
    ### Ordered variables first
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
    length.uniq = function(x) {
      length(unique(x))
    }
    cons = apply(dir, 1, length.uniq)

    # consist= (cons==1 & pval.comp > 0.05)
    signsum = function(x) {
      sum(sign(x))
    }
    consist = cons == 1 & abs(apply(theta, 1, signsum)) == V
    procs <- c("Holm", "BH")
    if (n > 1) {
      res <- mt.rawp2adjp(pvalue, procs)
      oo <- res$index
      outres = data.frame(factor = factr[oo], theta[oo, ], psi[oo],
                          CI95[oo], res$adj, lbs[oo, ], consist[oo])
    }
    if (n == 1) {
      outres = data.frame(factor = factr, theta, psi, CI95, rawp = pvalue,
                          Holm = pvalue, BH = pvalue, lbs, consist)
    }

    outres = outres[is.na(outres[, "rawp"]) == F, ]
    names(outres)[1:(2 * V)] = c("VarType", paste("psiV", 1:V, sep = ""),
                                 "AvePsi", "CI95")
    names(outres)[(9 + 2 * V + 1)] = "Consistent"

    ### Get Consistency Measure and only significant
    ### Make BH cut-off flexible in future versions (default at 0.05)
    outres.cons = outres[outres[, "BH"] < 0.05, , drop = F]
    outres.cons = outres.cons[outres.cons[, "Consistent"] == TRUE,
                              c("VarType", "AvePsi", "rawp"), drop = F]


    # drops = c('VarType','description','Holm,')
    # outres.all=outres[,!(names(outres) %in% drops)]
    outres.byV = outres[, c(2:(2 + V - 1), 9:(9 + 2 * V)), drop = F]

    print(xtable(outres.byV, caption = "data Variable Importance Results By Estimation Sample",
                 label = "byV", digits = 4), type = "latex", file = paste(dirout,
                                                                          outname, "byV.tex", sep = ""), caption.placement = "top", include.rownames = T)

    outres.all = outres[, c("AvePsi", "CI95", "rawp", "BH"), drop = F]
    print(xtable(outres.all, caption = "data Variable Importance Results for Combined Estimates",
                 label = "allRes", digits = 4), type = "latex", file = paste(dirout,
                                                                             outname, "AllReslts.tex", sep = ""), caption.placement = "top",
          include.rownames = T)


    print(xtable(outres.cons[, ], caption = "Subset of of Significant and ``Consistent'' Results",
                 label = "consisRes", digits = 4), type = "latex", file = paste(dirout,
                                                                                outname, "ConsistReslts.tex", sep = ""), caption.placement = "top",
          include.rownames = T)
  }

  ###########################################################################
  # end of just continuous variables
  ###########################################################################

  if (n.num > 0 & n.fac > 0) {
    ## For each factor, apply qq.range function and get rid of those where
    ## 'true' data.fac is data frame of variables that are factors
    xp = dim(data.fac)[2]
    facnames = names(data.fac)
    nam.fac = function(x, name) {
      nc = nchar(x)
      out = paste(name, substr(x, 2, nc), sep = "XX")
      out = gsub(" ", "", out)
      return(out)
    }
    newX = NULL
    coln = NULL
    options(na.action = "na.pass")
    for (i in 1:xp) {
      # cat(' i = ',i,'\n')
      x = data.fac[, i]
      inds <- model.matrix(~x - 1)[, -1]
      nmes = colnames(inds)
      if (is.null(nmes)) {
        nmes2 = facnames[i]
      }
      if (!is.null(nmes)) {
        nmes2 = nam.fac(nmes, facnames[i])
      }
      coln = c(coln, nmes2)
      newX = cbind(newX, inds)
    }
    colnames(newX) = coln
    ## Indexing vector for dummy basis back to original factors
    cc = regexpr("XX", coln)
    ncc = nchar(coln)
    cc[cc < 0] = ncc[cc < 0] + 1
    fac.indx = substr(coln, 1, cc - 1)
    datafac.dum = newX
    ### Now, make deciles for continuous variables
    X = data.num
    xc = dim(X)[2]
    qt = apply(na.omit(X), 2, quantile, probs = seq(0.1, 0.9, 0.1))
    newX = NULL
    coln = NULL
    varn = colnames(X)
    num.cat = apply(X, 2, ln.unique)
    ### Processing continuous variables
    Xnew = NULL
    for (k in 1:xc) {
      Xt = X[, k]
      tst = as.numeric(discretize(Xt, method = "frequency", categories = 10,
                                  ordered = T))
      Xnew = cbind(Xnew, tst)
    }
    colnames(Xnew) = varn
    data.cont.dist = Xnew
    ############### Missing Basis
    xp = dim(data.cont.dist)[2]
    n.cont = dim(data.cont.dist)[1]
    sum.na = function(x) {
      sum(is.na(x))
    }
    sna = apply(data.cont.dist, 2, sum.na)
    nmesX = colnames(data.cont.dist)
    miss.cont = NULL
    nmesm = NULL
    for (k in 1:xp) {
      if (sna[k] > 0) {
        ix = as.numeric(is.na(data.cont.dist[, k]) == F)
        miss.cont = cbind(miss.cont, ix)
        nmesm = c(nmesm, paste("Imiss_", nmesX[k], sep = ""))
      }
    }
    # if(is.null(miss.cont)){miss.cont= rep(1,n.cont)}
    colnames(miss.cont) = nmesm
    ############ Missing Basis for Factors
    xp = dim(datafac.dum)[2]
    sna = apply(datafac.dum, 2, sum.na)
    nmesX = colnames(datafac.dum)
    miss.fac = NULL
    nmesm = NULL
    for (k in 1:xp) {
      if (sna[k] > 0) {
        ix = as.numeric(is.na(datafac.dum[, k]) == F)
        datafac.dum[is.na(datafac.dum[, k]), k] = 0
        miss.fac = cbind(miss.fac, ix)
        nmesm = c(nmesm, paste("Imiss_", nmesX[k], sep = ""))
      }
    }
    colnames(miss.fac) = nmesm

    #### Start Estimator First using All Data Start with continuous variables

    numcat.cont = apply(data.cont.dist, 2, ln.unique)
    xc = length(numcat.cont)
    cats.cont <- lapply(1:xc, function(i) {
      sort(unique(data.cont.dist[, i]))
    })


    ## Find the level of covariate that has lowest risk
    datafac.dumW = datafac.dum
    datafac.dumW[is.na(datafac.dum)] = 0
    data.numW = data.num
    data.numW[is.na(data.num)] = 0
    n = length(Y)
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
      tmle.1 = tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                    Q.SL.library = Q.lib, family = fam, verbose = FALSE)
      g1 = tmle.1$g$g1W
      Qst = tmle.1$Qstar[, 2]
      theta = mean(Qst)
      IC = (A/g1) * (Y - Qst) + Qst - theta
      return(list(theta = theta, IC = IC))
    }

    #### Stratified CV to insure balance (by one grouping variable, Y)
    CC.CV = function(V, Y) {
      nn = length(Y)
      tt = table(Y)
      Ys = unique(Y)
      nys = length(Ys)
      out = rep(NA, nn)
      for (i in 1:nys) {
        n = as.numeric(tt[i])
        xx = cvFolds(n, K = V, R = 1, type = "random")$which
        out[Y == Ys[i]] = xx
      }
      return(out)
    }
    folds = CC.CV(V, Y)
    max.2 = function(x) {
      max(x^2, na.rm = T)
    }

    # detach('package:cvTools', unload=TRUE) detach('package:lattice',
    # unload=TRUE) library(doParallel) cl <- makeCluster(10)
    # registerDoParallel(cl) Below is to get indexing vectors so that any
    # basis functions related to current A that are in covariate matrix can
    # be removed
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
    xc = dim(data.fac)[2]
    n.fac = dim(data.fac)[1]
    ############################################################################### VIM for Factors
    output <- lapply(1:xc, function(i) {
      # output <- lapply(1:2, function(i) { pp1=proc.time()
      nameA = names.fac[i]
      cat(" i = ", i, "Var = ", nameA, " out of ", xc, " factor variables",
          "\n")
      thetaV = NULL
      varICV = NULL
      labV = NULL
      EY0V = NULL
      EY1V = NULL
      nV = NULL
      for (kk in 1:V) {
        cat(" i = ", i, " V = ", kk, "\n")
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
        if (is.null(missdumW)) {
          missdumW = rep(NA, n.fac)
        }
        if (is.null(miss.cont)) {
          miss.cont = rep(NA, n.fac)
        }
        if (is.null(dumW)) {
          dumW = rep(NA, n.fac)
        }
        if (is.null(data.numW)) {
          data.numW = rep(NA, n.fac)
        }
        W = data.frame(data.numW, miss.cont, dumW, missdumW)
        W = W[, !apply(is.na(W), 2, all), drop = F]
        Wt = W[folds != kk, , drop = F]
        Wv = W[folds == kk, , drop = F]
        Adum = data.frame(Adum[folds != kk, ])
        ### Pull out any variables that are overly correlated with At (corr coef
        ### < corthes)
        corAt = apply(cor(Adum, Wt, use = "complete.obs"), 2, max.2)
        corAt[corAt < -1] = 0
        # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
        incc = abs(corAt) < corthres & is.na(corAt) == F
        Wv = Wv[, incc, drop = F]
        Wt = Wt[, incc, drop = F]
        nw = dim(Wv)[2]
        ### Skip of number covariates < 10
        if (nw <= 10) {
          Wtsht = Wt
          Wvsht = Wv
        }
        if (nw > 10) {
          mydist <- as.matrix(distancematrix(t(Wt), d = "cosangle",
                                             na.rm = T))
          hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "mean",
                                 verbose = FALSE, K = 10, kmax = 3, khigh = 3), silent = TRUE)
          if (class(hopach.1) == "try-error") {
            hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "med",
                                   verbose = FALSE, K = 10, kmax = 3, khigh = 3), silent = TRUE)
          }
          if (class(hopach.1) == "try-error") {
            Wtsht = Wt
            Wvsht = Wv
          }
          if (class(hopach.1) != "try-error") {
            nlvls = nchar(max(hopach.1$final$labels))
            no <- trunc(mean(log10(hopach.1$final$labels)))
            ### Find highest level of tree where minimum number of covariates is >
            ### ncov
            lvl = 1:nlvls
            ncv = NULL
            for (ii in lvl) {
              ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no -
                                                                           (ii - 1))))))
            }
            ncv = unique(ncv)
            lev = min(min(nlvls, dim(Wt)[2]), min(lvl[ncv >= ncov]))
            two.clust <- unique(trunc(hopach.1$final$labels/(10^(no -
                                                                   (lev - 1)))))
            md <- hopach.1$final$medoids
            mm = md[, 1] %in% two.clust
            incc = md[mm, 2]
            Wtsht = Wt[, incc]
            Wvsht = Wv[, incc]
          }
        }
        # cat(' time a = ',proc.time()-pp1,'\n') pp2=proc.time()
        deltat = as.numeric(is.na(Yt) == F & is.na(At) == F)
        deltav = as.numeric(is.na(Yv) == F & is.na(Av) == F)
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
        vals = levA[pmin(minc, minc2) > minCell]
        num.cat = length(vals)
        nYt = sum(Yt[is.na(At) == FALSE])
        nYv = sum(Yv[is.na(Av) == FALSE])
        ## Don't do if 1) no more than one category of A left or if missingness
        ## pattern for A is such that there are few death events left in either
        ## (< minYs)
        if (num.cat < 2 | min(nYt, nYv) < minYs) {
          thetaV = c(thetaV, NA)
          varICV = c(varICV, NA)
          labV = rbind(labV, c(NA, NA))
          EY0V = c(EY0V, NA)
          EY1V = c(EY1V, NA)
          nV = c(nV, NA)
        }
        if (num.cat >= 2 & min(nYt, nYv) >= minYs) {
          labmin = NULL
          labmax = NULL
          maxEY1 = -1e+05
          minEY1 = 1e+06
          minj = 0
          maxj = 0
          ICmin = NULL
          ICmax = NULL
          errcnt = 0
          for (j in 1:num.cat) {
            # cat(' j = ',j,'\n')
            IA = as.numeric(At == vals[j])
            IA[is.na(IA)] = 0
            # if(min(table(IA,Yt))>=)
            res = try(get.tmle.est(Yt, IA, Wtsht, deltat, Q.lib = Q.library,
                                   g.lib = g.library), silent = TRUE)
            if (class(res) == "try-error") {
              errcnt = errcnt + 1
            }
            if (class(res) != "try-error") {
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
    })

    ###############################################################################
    # Repeat for Continuous Explanatory Variables
    ###############################################################################

    cor.two = function(x, y) {
      (cor(na.omit(cbind(x, y)))[1, 2])^2
    }
    folds = CC.CV(V, Y)
    # detach('package:cvTools', unload=TRUE) detach('package:lattice',
    # unload=TRUE) library(doParallel) cl <- makeCluster(20)
    # registerDoParallel(cl)
    names.cont = colnames(data.cont.dist)
    xc = dim(data.cont.dist)[2]
    n.cont = dim(data.cont.dist)[1]
    ######################################################
    out.put <- lapply(1:xc, function(i) {
      nameA = names.cont[i]
      cat(" i = ", i, "Var = ", nameA, " out of ", xc, " numeric variables",
          "\n")
      thetaV = NULL
      varICV = NULL
      labV = NULL
      EY0V = NULL
      EY1V = NULL
      nV = NULL
      for (kk in 1:V) {
        # cat(' v = ',kk,' out of V = ',V,'\n')
        At = data.cont.dist[folds != kk, i]
        Av = data.cont.dist[folds == kk, i]
        Yt = Y[folds != kk]
        Yv = Y[folds == kk]
        AY1 = At[Yt == 1 & is.na(At) == F]
        hh = histogram::histogram(AY1, verbose = F, type = "irregular")$breaks
        if (hh[length(hh)] < max(At, na.rm = T)) {
          hh[length(hh)] = max(At, na.rm = T) + 0.1
        }
        if (hh[1] > min(At[At > 0], na.rm = T)) {
          hh[1] = min(At[At > 0], na.rm = T) - 0.1
        }
        Atnew = cut(At, breaks = hh)
        Avnew = cut(Av, breaks = hh)
        if (length(na.omit(unique(Atnew))) <= 1 | length(na.omit(unique(Avnew))) <=
            1) {
          thetaV = c(thetaV, NA)
          varICV = c(varICV, NA)
          labV = rbind(labV, c(NA, NA))
          EY0V = c(EY0V, NA)
          EY1V = c(EY1V, NA)
          nV = c(nV, NA)
        }
        if (length(na.omit(unique(Atnew))) > 1 & length(na.omit(unique(Avnew))) >
            1) {
          labs = names(table(Atnew))
          Atnew = as.numeric(Atnew) - 1
          Avnew = as.numeric(Avnew) - 1
          numcat.cont[i] = length(labs)
          # change this to match what was done for factors - once
          # cats.cont[[i]]=as.numeric(na.omit(unique(Atnew)))
          cats.cont[[i]] = as.numeric(names(table(Atnew)))
          ### acit.numW is just same as data.cont.dist except with NA's replaced by
          ### 0's.
          if (is.null(miss.cont)) {
            miss.cont = rep(NA, n.cont)
          }
          if (is.null(miss.fac)) {
            miss.fac = rep(NA, n.cont)
          }
          if (is.null(datafac.dumW)) {
            datafac.dumW = rep(NA, n.cont)
          }
          if (is.null(data.numW)) {
            data.numW = rep(NA, n.cont)
          }
          W = data.frame(data.numW[, -i, drop = F], miss.cont,
                         datafac.dumW, miss.fac)
          W = W[, !apply(is.na(W), 2, all), drop = F]
          Wt = W[folds != kk, , drop = F]
          Wv = W[folds == kk, , drop = F]

          nmesW = names(Wt)
          mtch = match(nmesW, paste("Imiss_", nameA, sep = ""))
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
          nw = dim(Wv)[2]
          ### Skip of number covariates < 10
          if (nw <= 10) {
            Wtsht = Wt
            Wvsht = Wv
          }
          if (nw > 10) {
            mydist <- as.matrix(distancematrix(t(Wt), d = "cosangle",
                                               na.rm = T))
            hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "mean",
                                   verbose = FALSE), silent = TRUE)
            if (class(hopach.1) == "try-error") {
              hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "med",
                                     verbose = FALSE), silent = TRUE)
            }
            if (class(hopach.1) == "try-error") {
              Wtsht = Wt
              Wvsht = Wv
            }
            if (class(hopach.1) != "try-error") {
              nlvls = nchar(max(hopach.1$final$labels))
              no <- trunc(mean(log10(hopach.1$final$labels)))
              ### Find highest level of tree where minimum number of covariates is >
              ### ncov
              lvl = 1:nlvls
              ncv = NULL
              for (ii in lvl) {
                ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no -
                                                                             (ii - 1))))))
              }
              ncv = unique(ncv)
              lev = min(min(nlvls, dim(Wt)[2]), min(lvl[ncv >=
                                                          ncov]))
              two.clust <- unique(trunc(hopach.1$final$labels/(10^(no -
                                                                     (lev - 1)))))
              md <- hopach.1$final$medoids
              mm = md[, 1] %in% two.clust
              incc = md[mm, 2]
              Wtsht = Wt[, incc]
              Wvsht = Wv[, incc]
            }
          }
          deltat = as.numeric(is.na(Yt) == F & is.na(Atnew) ==
                                F)
          deltav = as.numeric(is.na(Yv) == F & is.na(Avnew) ==
                                F)
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
          if (min(table(Avnew[Avnew >= 0], Yv[Avnew >= 0])) <=
              minCell) {
            thetaV = c(thetaV, NA)
            varICV = c(varICV, NA)
            labV = rbind(labV, c(NA, NA))
            EY0V = c(EY0V, NA)
            EY1V = c(EY1V, NA)
            nV = c(nV, NA)
          }
          if (min(table(Avnew, Yv)) > minCell) {
            labmin = NULL
            labmax = NULL
            errcnt = 0
            for (j in 1:numcat.cont[i]) {
              # cat(' i = ',i,' kk = ',kk,' j = ',j,'\n')
              IA = as.numeric(Atnew == vals[j])
              res = try(get.tmle.est(Yt, IA, Wtsht, deltat, Q.lib = Q.library,
                                     g.lib = g.library), silent = TRUE)
              if (class(res) == "try-error") {
                errcnt = errcnt + 1
              }
              if (class(res) != "try-error") {
                EY1 = res$theta
                if (EY1 < minEY1) {
                  minj = j
                  minEY1 = EY1
                  labmin = labs[j]
                }
                if (EY1 > maxEY1) {
                  maxj = j
                  maxEY1 = EY1
                  labmax = labs[j]
                }
              }
            }
            ##### Now, estimate on validation sample
            if (errcnt == numcat.cont[i] | minj == maxj) {
              thetaV = c(thetaV, NA)
              varICV = c(varICV, NA)
              labV = rbind(labV, c(NA, NA))
              EY0V = c(EY0V, NA)
              EY1V = c(EY1V, NA)
              nV = c(nV, NA)
            }
            if (errcnt != numcat.cont[i] & minj != maxj) {
              IA = as.numeric(Avnew == vals[minj])
              res = try(get.tmle.est(Yv, IA, Wvsht, deltav, Q.lib = Q.library,
                                     g.lib = g.library), silent = TRUE)
              if (class(res) == "try-error") {
                thetaV = c(thetaV, NA)
                varICV = c(varICV, NA)
                labV = rbind(labV, c(NA, NA))
                EY0V = c(EY0V, NA)
                EY1V = c(EY1V, NA)
                nV = c(nV, NA)
              }
              if (class(res) != "try-error") {
                IC0 = res$IC
                EY0 = res$theta
                IA = as.numeric(Avnew == vals[maxj])
                res2 = try(get.tmle.est(Yv, IA, Wvsht, deltav,
                                        Q.lib = Q.library, g.lib = g.library), silent = TRUE)
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
      }
      list(EY1V, EY0V, thetaV, varICV, labV, nV, "numeric")
    })

    vars.cont = colnames(data.cont.dist)
    out.cont = out.put
    stf = c("vars.cont", "out.cont")
    # allt=ls() mm=allt%in%stf allt=allt[mm==F] rm(list=allt)
    vars = colnames(data.fac)
    nc = length(vars.cont)
    nf = length(vars)
    factr = c(rep("ordered", nc), rep("factor", nf))
    vars = c(vars.cont, vars)
    out.put = c(out.cont, output)
    names(out.put) = vars
    lngth = sapply(out.put, function(x) length(x))
    out.put = out.put[lngth == 7]
    vars = vars[lngth == 7]
    factr = factr[lngth == 7]
    ## Get rid of any variables that have a validation sample with no
    ## estimates of variable importance
    if (length(out.put) == 0) {
      write("No VIM's could be calculated due to sample size, etc",
            file = "AllReslts.tex")
    }
    if (length(out.put) > 0) {
      lngth2 = sapply(out.put, function(x) length(na.omit(x[[1]])))
      out.put = out.put[lngth2 == V]
      vars = vars[lngth2 == V]
      factr = factr[lngth2 == V]
      if (length(out.put) == 0) {
        write("No VIM's could be calculated due to sample size, etc",
              file = "AllReslts.tex")
      }
      if (length(out.put) > 0) {
        tst = lapply(out.put, function(x) x[[3]])
        tst = do.call(rbind, tst)
        tot.na = function(x) {
          sum(is.na(x))
        }
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
          names(out) = paste("v.", lbel, rep(c("a_L", "a_H"), 2),
                             sep = "")
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
        CI95 = paste("(", signif(lower, signdig), " - ", signif(upper,
                                                                signdig), ")", sep = "")
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
        consist = cons == 1 & abs(apply(theta, 1, signsum)) ==
          V
        procs <- c("Holm", "BH")
        if (n > 1) {
          res <- mt.rawp2adjp(pvalue, procs)
          oo <- res$index
          outres = data.frame(factor = factr[oo], theta[oo, ],
                              psi[oo], CI95[oo], res$adj, lbs[oo, ], consist[oo])
        }
        if (n == 1) {
          outres = data.frame(factor = factr, theta, psi, CI95,
                              rawp = pvalue, Holm = pvalue, BH = pvalue, lbs, consist)
        }
        outres = outres[is.na(outres[, "rawp"]) == F, , drop = F]
        names(outres)[1:(1 + 2 * V)] = c("VarType", paste("psiV",
                                                          1:V, sep = ""), "AvePsi", "CI95")
        names(outres)[(9 + 2 * V)] = "Consistent"

        ### Get Consistency Measure and only significant Make BH cut-off flexible
        ### in future versions (default at 0.05)
        outres.cons = outres[outres[, "BH"] < 0.05, , drop = F]
        outres.cons = outres.cons[outres.cons[, "Consistent"] ==
                                    TRUE, c("VarType", "AvePsi", "rawp"), drop = F]


        # drops = c('VarType','description','Holm,')
        # outres.all=outres[,!(names(outres) %in% drops)]
        outres.byV = outres[, c(2:(2 + V - 1), 9:(9 + 2 * V)),
                            drop = F]
        outres.all = outres[, c("AvePsi", "CI95", "rawp", "BH"),
                            drop = F]

        print(xtable(outres.byV, caption = "data Variable Importance Results By Estimation Sample",
                     label = "byV", digits = 4), type = "latex", file = paste(dirout,
                                                                              outname, "byV.tex", sep = ""), caption.placement = "top",
              include.rownames = T)

        print(xtable(outres.all, caption = "data Variable Importance Results for Combined Estimates",
                     label = "allRes", digits = 4), type = "latex", file = paste(dirout,
                                                                                 outname, "AllReslts.tex", sep = ""), caption.placement = "top",
              include.rownames = T)


        print(xtable(outres.cons, caption = "Subset of of Significant and ``Consistent'' Results",
                     label = "consisRes", digits = 4), type = "latex", file = paste(dirout,
                                                                                    outname, "ConsistReslts.tex", sep = ""), caption.placement = "top",
              include.rownames = T)
      }
    }
  }

  ###############################################################
  # If only have factors for predictors
  ###############################################################

  if (n.fac > 0 & n.num == 0) {
    ## For each factor, apply qq.range function and get rid of those where
    ## 'true' data.fac is data frame of variables that are factors
    xp = dim(data.fac)[2]
    facnames = names(data.fac)
    nam.fac = function(x, name) {
      nc = nchar(x)
      out = paste(name, substr(x, 2, nc), sep = "XX")
      out = gsub(" ", "", out)
      return(out)
    }
    newX = NULL
    coln = NULL
    options(na.action = "na.pass")
    for (i in 1:xp) {
      # cat(' i = ',i,'\n')
      x = data.fac[, i]
      inds <- model.matrix(~x - 1)[, -1]
      nmes = colnames(inds)
      if (is.null(nmes)) {
        nmes2 = facnames[i]
      }
      if (!is.null(nmes)) {
        nmes2 = nam.fac(nmes, facnames[i])
      }
      coln = c(coln, nmes2)
      newX = cbind(newX, inds)
    }
    colnames(newX) = coln
    ## Indexing vector for dummy basis back to original factors
    cc = regexpr("XX", coln)
    ncc = nchar(coln)
    cc[cc < 0] = ncc[cc < 0] + 1
    fac.indx = substr(coln, 1, cc - 1)
    datafac.dum = newX
    ############ Missing Basis for Factors
    xp = dim(datafac.dum)[2]
    sna = apply(datafac.dum, 2, sum.na)
    nmesX = colnames(datafac.dum)
    miss.fac = NULL
    nmesm = NULL
    for (k in 1:xp) {
      if (sna[k] > 0) {
        ix = as.numeric(is.na(datafac.dum[, k]) == F)
        datafac.dum[is.na(datafac.dum[, k]), k] = 0
        miss.fac = cbind(miss.fac, ix)
        nmesm = c(nmesm, paste("Imiss_", nmesX[k], sep = ""))
      }
    }
    colnames(miss.fac) = nmesm

    ## Find the level of covariate that has lowest risk
    datafac.dumW = datafac.dum
    datafac.dumW[is.na(datafac.dum)] = 0
    n = length(Y)
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
      tmle.1 = tmle(Y, A, W, Delta = delta, g.SL.library = g.lib,
                    Q.SL.library = Q.lib, family = fam, verbose = FALSE)
      g1 = tmle.1$g$g1W
      Qst = tmle.1$Qstar[, 2]
      theta = mean(Qst)
      IC = (A/g1) * (Y - Qst) + Qst - theta
      return(list(theta = theta, IC = IC))
    }

    #### Stratified CV to insure balance (by one grouping variable, Y)
    CC.CV = function(V, Y) {
      nn = length(Y)
      tt = table(Y)
      Ys = unique(Y)
      nys = length(Ys)
      out = rep(NA, nn)
      for (i in 1:nys) {
        n = as.numeric(tt[i])
        xx = cvFolds(n, K = V, R = 1, type = "random")$which
        out[Y == Ys[i]] = xx
      }
      return(out)
    }
    folds = CC.CV(V, Y)
    max.2 = function(x) {
      max(x^2, na.rm = T)
    }

    # detach('package:cvTools', unload=TRUE) detach('package:lattice',
    # unload=TRUE) library(doParallel) cl <- makeCluster(10)
    # registerDoParallel(cl) Below is to get indexing vectors so that any
    # basis functions related to current A that are in covariate matrix can
    # be removed
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
    xc = dim(data.fac)[2]
    n.fac = dim(data.fac)[1]

    ###############################################################################
    # VIM for Factors
    ###############################################################################

    output <- lapply(1:xc, function(i) {
      nameA = names.fac[i]
      cat(" i = ", i, "Var = ", nameA, " out of ", xc, " factor variables",
          "\n")
      thetaV = NULL
      varICV = NULL
      labV = NULL
      EY0V = NULL
      EY1V = NULL
      nV = NULL
      for (kk in 1:V) {
        # cat(' i = ',i,' V = ',kk,'\n')
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
        if (is.null(missdumW)) {
          missdumW = rep(NA, n.fac)
        }
        if (is.null(dumW)) {
          dumW = rep(NA, n.fac)
        }
        W = data.frame(dumW, missdumW)
        W = W[, !apply(is.na(W), 2, all), drop = F]
        Wt = W[folds != kk, , drop = F]
        Wv = W[folds == kk, , drop = F]
        Adum = data.frame(Adum[folds != kk, ])
        ### Pull out any variables that are overly correlated with At (corr coef
        ### < corthes)
        corAt = apply(cor(Adum, Wt, use = "complete.obs"), 2, max.2)
        # cat('i = ',i,' maxCor = ',max(corAt,na.rm=T),'\n')
        incc = abs(corAt) < corthres & is.na(corAt) == F
        Wv = Wv[, incc, drop = F]
        Wt = Wt[, incc, drop = F]
        nw = dim(Wv)[2]
        ### Skip of number covariates < 10
        if (nw <= 10) {
          Wtsht = Wt
          Wvsht = Wv
        }
        if (nw > 10) {
          mydist <- as.matrix(distancematrix(t(Wt), d = "cosangle",
                                             na.rm = T))
          hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "mean",
                                 verbose = FALSE), silent = TRUE)
          if (class(hopach.1) == "try-error") {
            hopach.1 <- try(hopach(t(Wt), dmat = mydist, mss = "med",
                                   verbose = FALSE), silent = TRUE)
          }
          if (class(hopach.1) == "try-error") {
            Wtsht = Wt
            Wvsht = Wv
          }
          if (class(hopach.1) != "try-error") {
            nlvls = nchar(max(hopach.1$final$labels))
            no <- trunc(mean(log10(hopach.1$final$labels)))
            ### Find highest level of tree where minimum number of covariates is >
            ### ncov
            lvl = 1:nlvls
            ncv = NULL
            for (ii in lvl) {
              ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no -
                                                                           (ii - 1))))))
            }
            ncv = unique(ncv)
            lev = min(min(nlvls, dim(Wt)[2]), min(lvl[ncv >= ncov]))
            two.clust <- unique(trunc(hopach.1$final$labels/(10^(no -
                                                                   (lev - 1)))))
            md <- hopach.1$final$medoids
            mm = md[, 1] %in% two.clust
            incc = md[mm, 2]
            Wtsht = Wt[, incc]
            Wvsht = Wv[, incc]
          }
        }

        deltat = as.numeric(is.na(Yt) == F & is.na(At) == F)
        deltav = as.numeric(is.na(Yv) == F & is.na(Av) == F)
        if (sum(deltat == 0) < 10) {
          Yt = Yt[deltat == 1]
          At = At[deltat == 1]
          Wtsht = Wtsht[deltat == 1, ]
          deltat = deltat[deltat == 1]
        }
        levA = levels(At)
        minc = apply(table(Av, Yv), 1, min)
        minc2 = apply(table(At, Yt), 1, min)
        vals = levA[pmin(minc, minc2) > minCell]
        num.cat = length(vals)
        nYt = sum(Yt[is.na(At) == FALSE])
        nYv = sum(Yv[is.na(Av) == F])
        ## Don't do if 1) no more than one category of A left or if missingness
        ## pattern for A is such that there are few death events left in either
        ## (< minYs)
        if (num.cat < 2 | min(nYt, nYv) < minYs) {
          thetaV = c(thetaV, NA)
          varICV = c(varICV, NA)
          labV = rbind(labV, c(NA, NA))
          EY0V = c(EY0V, NA)
          EY1V = c(EY1V, NA)
          nV = c(nV, NA)
        }
        if (num.cat >= 2 & min(nYt, nYv) >= minYs) {
          labmin = NULL
          labmax = NULL
          maxEY1 = -1e+05
          minEY1 = 1e+06
          minj = 0
          maxj = 0
          ICmin = NULL
          ICmax = NULL
          errcnt = 0
          for (j in 1:num.cat) {
            cat(" j = ", j, "\n")
            IA = as.numeric(At == vals[j])
            IA[is.na(IA)] = 0
            # if(min(table(IA,Yt))>=)
            res = try(get.tmle.est(Yt, IA, Wtsht, deltat, Q.lib = Q.library,
                                   g.lib = g.library))
            if (class(res) == "try-error") {
              errcnt = errcnt + 1
            }
            if (class(res) != "try-error") {
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
          ##### Now, estimate on validation sample
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
                                   g.lib = g.library))
            if (class(res) == "try-error") {
              thetaV = c(thetaV, NA)
              varICV = c(varICV, NA)
              labV = rbind(labV, c(NA, NA))
              EY0V = c(EY0V, NA)
              EY1V = c(EY1V, NA)
              nV = c(nV, NA)
            }
            if (class(res) != "try-error") {
              IC0 = res$IC
              EY0 = res$theta
              IA = as.numeric(Av == vals[maxj])
              IA[is.na(IA)] = 0
              res2 = try(get.tmle.est(Yv, IA, Wvsht, deltav, Q.lib = Q.library,
                                      g.lib = g.library))
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
      list(EY1V, EY0V, thetaV, varICV, labV, nV)
      # print(data.frame(EY1V,EY0V,thetaV,varICV,labV,nV))
    })

    ###############################################################################
    # Repeat for Continuous Explanatory Variables
    ###############################################################################

    vars = colnames(data.fac)
    # nc=length(vars.cont)
    nf = length(vars)
    factr = rep("factor", nf)
    out.put = output
    names(out.put) = vars
    lngth = sapply(out.put, function(x) length(x))
    out.put = out.put[lngth == 6]
    vars = vars[lngth == 6]
    factr = factr[lngth == 6]
    lngth2 = sapply(out.put, function(x) length(na.omit(x[[1]])))
    out.put = out.put[lngth2 == V]
    vars = vars[lngth2 == V]
    factr = factr[lngth2 == V]

    tst = lapply(out.put, function(x) x[[3]])
    tst = do.call(rbind, tst)
    tot.na = function(x) {
      sum(is.na(x))
    }
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
    tst = lapply(tst, labs.get, fold = 2)
    lbs = do.call(rbind, tst)


    meanvarIC = apply(varIC, 1, mean)
    psi = apply(theta, 1, mean)
    SE = sqrt(meanvarIC/n)
    lower = psi - 1.96 * SE
    upper = psi + 1.96 * SE
    signdig = 3

    # TODO: convert to paste0?
    CI95 = paste("(", signif(lower, signdig), " - ", signif(upper, signdig), ")", sep = "")

    # 1-sided p-value
    pvalue = 1 - pnorm(psi/SE)

    #####
    # FOR THETA (generalize to chi-square test?)
    # TT = (theta[,1]-theta[,2])/sqrt(SEV[,1]^2+SEV[,2]^2)
    # pval.comp=2*(1-pnorm(abs(TT))) FOR levels (just make sure in same order)
    nc = sum(factr == "ordered")
    n = length(factr)
    #### Factors
    lwr = NULL
    uwr = NULL
    for (i in 1:V) {
      lwr = cbind(lwr, lbs[(nc + 1):n, i * 2 - 1])
      uwr = cbind(uwr, lbs[(nc + 1):n, i * 2 - 1])
    }
    length.uniq = function(x) {
      length(unique(x))
    }
    conslwr = apply(lwr, 1, length.uniq)
    consupr = apply(uwr, 1, length.uniq)
    cons = conslwr * consupr
    # consist= (cons==1 & pval.comp < 0.05)
    signsum = function(x) {
      sum(sign(x))
    }
    consist = cons == 1 & abs(apply(theta, 1, signsum)) == V

    procs <- c("Holm", "BH")
    res <- mt.rawp2adjp(pvalue, procs)
    oo <- res$index
    outres = data.frame(factor = factr[oo], theta[oo, ], psi[oo], CI95[oo],
                        res$adj, lbs[oo, ], consist[oo])

    outres = outres[is.na(outres[, "rawp"]) == F, ]
    names(outres)[1:(1 + 2 * V)] = c("VarType", paste("psiV", 1:V,
                                                      sep = ""), "AvePsi", "CI95")
    names(outres)[(9 + 2 * V)] = "Consistent"

    #######
    # Get Consistency Measure and only significant
    # Make BH cut-off flexible in future versions (default at 0.05)
    outres.cons = outres[outres[, "BH"] < 0.05, , drop = F]
    outres.cons = outres.cons[outres.cons[, "Consistent"] == TRUE,
                              c("VarType", "AvePsi", "rawp"), drop = F]


    # drops = c('VarType','description','Holm,')
    # outres.all=outres[,!(names(outres) %in% drops)]
    outres.byV = outres[, c(2:(2 + V - 1), 9:(9 + 2 * V)), drop = F]

    print(xtable(outres.byV, caption = "data Variable Importance Results By Estimation Sample",
                 label = "byV", digits = 4), type = "latex", file = paste(dirout,
                                                                          outname, "byV.tex", sep = ""), caption.placement = "top", include.rownames = T)

    outres.all = outres[, c("AvePsi", "CI95", "rawp", "BH"), drop = F]
    print(xtable(outres.all, caption = "data Variable Importance Results for Combined Estimates",
                 label = "allRes", digits = 4), type = "latex", file = paste(dirout,
                                                                             outname, "AllReslts.tex", sep = ""), caption.placement = "top",
          include.rownames = T)


    print(xtable(outres.cons[, ], caption = "Subset of of Significant and ``Consistent'' Results",
                 label = "consisRes", digits = 4), type = "latex", file = paste(dirout,
                                                                                outname, "ConsistReslts.tex", sep = ""), caption.placement = "top",
          include.rownames = T)
  }
}


