vim_factors =
  function(Y,
           numerics,
           factors,
           V,
           folds,
           A_names,
           family,
           minCell,
           minYs,
           Q.library,
           g.library,
           Qbounds,
           corthres,
           adjust_cutoff,
           verbose = FALSE,
           verbose_tmle = FALSE,
           verbose_reduction = FALSE) {

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

    xc = ncol(factors$data.fac)
    n.fac = nrow(factors$data.fac)

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
        At = factors$data.fac[folds != fold_k, var_i]

        # All data in this fold is the validation data.
        Av = factors$data.fac[folds == fold_k, var_i]

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
          Wtsht = Wtsht[deltat == 1, , drop = FALSE]
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
        if ((length(unique(Yt)) == 2L &&
             (num.cat < 2L || min(nYt, nYv) < minYs)) ||
            (length(is_constant) > 0 && mean(is_constant) == 1)) {
          if (length(is_constant) > 0 && mean(is_constant) == 1) {
            error_msg = paste("Skipping", nameA, "because HOPACH reduced W to",
                              "all constant columns.")
          } else if (num.cat < 2L) {
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

    colnames_factor = colnames(factors$data.fac)
  } else {
    colnames_factor = NULL
    vim_factor = NULL
    cat("No factor variables - skip VIM estimation.\n\n")
  }

  # Compile and return results.
  (results = list(
    vim_factor = vim_factor,
    colnames_factor = colnames_factor
  ))
}
