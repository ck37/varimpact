vim_numerics =
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
  # TODO: move out.
  cor.two = function(x, y) {
    (stats::cor(na.omit(cbind(x, y)))[1, 2])^2
  }

  # We use && so that the second check will be skipped when num_numeric == 0.
  if (numerics$num_numeric > 0 && ncol(numerics$data.cont.dist) > 0) {
    cat("Estimating variable importance for", numerics$num_numeric, "numerics.\n")

    xc = ncol(numerics$data.cont.dist)
    names.cont = colnames(numerics$data.cont.dist)
    n.cont = nrow(numerics$data.cont.dist)

    # Tally the number of unique values (bins) in each numeric variable; save as a vector.
    # ignore NAs though.
    numcat.cont = apply(numerics$data.cont.dist, 2, length_unique)

    if (verbose) {
      cat("Unique values by variable:\n")
      print(numcat.cont)
    }

    cats.cont = lapply(1:xc, function(i) {
      sort(unique(numerics$data.cont.dist[, i]))
    })

    ### Loop over each numeric variable.
    # Define var_i just to avoid automated NOTEs, will be overwritten by foreach.
    var_i = NULL
    #vim_numeric = foreach::foreach(var_i = 1:num_numeric, .verbose = verbose,
    #                               .errorhandling = "stop") %do_op% {
    vim_numeric = future.apply::future_lapply(1:numerics$num_numeric, future.seed = TRUE, function(var_i) {
    # TODO: reenable future_lapply
    #vim_numeric = lapply(1:numerics$num_numeric, function(var_i) {
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
          # To hold the level-specific results.
          levels = list(),
          # A dataframe to compile the level-specific results.
          bin_df = NULL,
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
        # Pr(A | Y = 1).
        if (length(unique(Yt)) == 2L) {
          # Binary outcome.

          # Conditional distribution of A given Y = 1.
          # F(A | Y = 1)
          A_Y1 = At[Yt == 1 & !is.na(At)]

          # Check if AY1 has only a single observation. If so, skip histogramming to avoid an error.
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
          penalized_hist = histogram::histogram(A_Y1, verbose = FALSE, type = "irregular", plot = FALSE)
          hist_breaks = penalized_hist$breaks

          if (verbose) {
            cat(nameA, "revised breaks:", hist_breaks, "\n")
          }

          # TODO: see if these next two steps are ever used/needed.

          # Check if the final cut-off is less that the maximum possible level; if so extend to slightly
          # larger than the maximimum possible level.
          if (hist_breaks[length(hist_breaks)] < max(At, na.rm = TRUE)) {
            hist_breaks[length(hist_breaks)] = max(At, na.rm = TRUE) + 0.1
          }

          # Check if the lowest cut-off is greater than the minimum possible bin; if so extend to slightly
          # below the minimum level.
          if (hist_breaks[1] > min(At[At > 0], na.rm = TRUE)) {
            hist_breaks[1] = min(At[At > 0], na.rm = TRUE) - 0.1
          }

          # Re-bin the training and validation vectors for the treatment variable based on the penalized
          # histogram.
          # This is creating factors, with levels specific to this CV-TMLE fold.
          Atnew = cut(At, breaks = hist_breaks)
          Avnew = cut(Av, breaks = hist_breaks)

          # TODO: check if the binning results in no-variation, and handle separately from the below situation.
          if (length(unique(Atnew)) == 1L) {
            warning(paste0(nameA, "has been converted to a single level after penalized histogramming."))
            #browser()
          }

        }
        if (singleAY1 || length(na.omit(unique(Atnew))) <= 1 ||
            length(na.omit(unique(Avnew))) <= 1) {
          error_msg = paste("Skipping", nameA, "in this fold because there is no variation.")
          if (verbose) cat(error_msg, "\n")
          fold_result$message = error_msg
          #warning(error_msg)
        } else {

          # These labels are simply the quantiles right now, not yet accounting
          # for the penalized histogramming.
          At_bin_labels = names(table(Atnew))

          # Non-discretized version of A in the training data; converted to a vector.
          At_raw = numerics$data.num[folds != fold_k, var_i]

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
          if (verbose) {
            cat("Updating numcat.cont for", nameA, "to be", length(At_bin_labels), "\n")
            print(At_bin_labels)
          }
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
            if (length(is_constant) > 0L) {
              cat("Constant columns (", length(is_constant), "):\n")
              print(is_constant)
            }
          }

          # Indicator that Y and A are both defined.
          deltat = as.numeric(!is.na(Yt) & !is.na(Atnew))
          deltav = as.numeric(!is.na(Yv) & !is.na(Avnew))

          # TODO: may want to remove this procedure, which is pretty arbitrary.
          if (sum(deltat == 0) < 10L) {
            Yt = Yt[deltat == 1]
            Wtsht = Wtsht[deltat == 1, , drop = FALSE]
            Atnew = Atnew[deltat == 1]
            deltat = deltat[deltat == 1]
          }

          vals = cats.cont[[var_i]]
          num.cat = length(vals)

          Atnew[is.na(Atnew)] = -1
          Avnew[is.na(Avnew)] = -1

          if ((length(is_constant) > 0 && mean(is_constant) == 1) ||
              (length(unique(Yt)) == 2L && min(table(Avnew[Avnew >= 0], Yv[Avnew >= 0])) <= minCell)) {
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

            # TODO: stop using training_estimates and switch to using bin_results.
            training_estimates = list()
            bin_results = list()

            ###########################################
            # Loop over each bin/level for this variable.

            # TODO: move into its own function?

            for (bin_j in 1:numcat.cont[var_i]) {

              # Create a list to hold the results for this level.
              bin_result = list(
                name = nameA,
                cv_fold = fold_k,
                level = bin_j,
                level_label = At_bin_labels[bin_j],

                # Training: Estimates
                # Our default value needs to be NA rather than NULL,
                # to allow rbinding into a dataframe later on.
                train_theta_tmle = NA,
                # TODO: implement these.
                #train_theta_iptw = NA,
                #train_theta_gcomp = NA,
                train_theta_unadj = NA,

                # Training: Misc
                train_cell_size = NA,
                train_msg = NA,

                # Test: Estimates
                test_theta_tmle = NA,
                # TODO: implement these.
                #test_theta_iptw = NA,
                #test_theta_gcomp = NA,
                test_theta_unadj = NA,

                # Test: Misc
                test_cell_size = NA,
                test_msg = NA,

                # Test: Predicted values (g, Q_bar, h)
                test_predictions = NULL
                #test_pred_g = NULL,
                #test_pred_Q_bar = NULL,
                #test_pred_h = NULL
              )


              # Create a treatment indicator, where 1 = obs in this bin
              # and 0 = obs not in this bin.
              IA = as.numeric(Atnew == vals[bin_j])

              # Save how many obs have this level/bin in this training fold.
              bin_result$train_cell_size = sum(IA)

              ###############################################
              # TODO: send to calculate_estimates first, which would calculate
              # TMLE, IPTW, G-Comp, and Unadj estimates.

              # Save unadjusted estimate: outcome mean among observations
              # at the desired treatment level, who are not missing their outcome value.
              bin_result$train_theta_unadj = mean(Yt[IA & deltat])

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
                tmle_result$label = At_bin_labels[bin_j]

                training_estimates[[bin_j]] = tmle_result

                # TODO: may also want to save the TMLE object to this bin_result list.
                bin_result$train_theta_tmle = tmle_result$theta

                if (verbose) cat(".")
              }

              #######################################################
              # NEW: also run code on corresponding validation fold.

              # TODO: remove the later validation code and use these results instead.

              # Indicator for having the desired treatment bin on validation
              IA = as.numeric(Avnew == vals[bin_j])

              # Missing values are not taken to be in this level.
              IA[is.na(IA)] = 0

              # Save how many obs have this level/bin in this validation fold.
              bin_result$test_cell_size = sum(IA)

              ##################
              # Run estimates on validation data (TMLE, IPTW, G-Comp, Unadj)
              # TODO: move into its own function.

              # Save unadjusted estimate: outcome mean among observations
              # at the desired treatment level, who are not missing their outcome value.
              bin_result$test_theta_unadj = mean(Yv[IA & deltav])

              # CV-TMLE: predict g, Q, and clever covariate on validation data.
              if (!is.null(training_estimates[[bin_j]])) {
                preds = try(apply_tmle_to_validation(Yv, IA, Wvsht, family,
                                                   deltav, training_estimates[[bin_j]],
                                                   verbose = verbose))
                if (class(preds) == "try-error") {
                  bin_result$test_msg = paste("CV-TMLE prediction on validation failed")
                } else {
                  # Save the result.
                  bin_result$test_predictions = preds
                }
              }

              bin_result$test_msg = "success"

              # Save to the main list.
              bin_results[[bin_j]] = bin_result
            }

            # Finished looping over each level of the assignment variable
            # (primarily training fold, but now also some val fold work).
            if (verbose) cat(" done.\n")

            # Save individual bin results.
            fold_result$bin_results = bin_results

            # Create a dataframe version of the bin results.
            fold_result$bin_df =
              do.call(rbind, lapply(bin_results, function(result) {
                # Exclude certain elements from the list - here the prediction vectors.
                # These should be saved separately.
                data.frame(result[!names(result) %in% c("test_predictions")],
                          stringsAsFactors = FALSE)
            }))

            # Save test_predictions for each bin into a combined dataframe.
            fold_result$test_predictions =
              do.call(rbind, lapply(1:length(bin_results), function(bin) {
                result = bin_results[[bin]]
                tryCatch({
                  data.frame(bin = bin,
                             bin_label = At_bin_labels[bin],
                             fold = fold_k,
                             result$test_predictions,
                             stringsAsFactors = FALSE)
                 }, error = function(error) {
                   NULL
                 })
              })
            )

            #####################################
            # Resume normal varimpact algorithm.

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
              # TODO: use the validation results already saved in bin_result

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

      ########################
      # Reconstruct CV-TMLE treatment-specific mean estimate for each validation fold.
      # Loop over bins and calculate bin-specific CV-TMLE estimate.
      # TODO: move this code into its own function
      if (T) {
      if (verbose) cat(paste0(nameA, ":"), "Estimating CV-TMLE treatment-specific means.\n")
      for (bin in 1:numcat.cont[var_i]) {

        if (verbose) {
          cat("Var:", var_i, "Bin", bin, "of", numcat.cont[var_i], "\n")
        }
        # Expects a list of results by fold.
        # Each element of that list should have the val_preds list, which
        # is calculated by apply_tmle_to_validation and currently saved in
        # fold_results[[*]]$test_predictions (separately by fold * level).
        combine_rows = lapply(fold_results, function(fold_r) {
          # Extract the rows specific to this bin/level.
          rows = fold_r$test_predictions[fold_r$test_predictions$bin == bin, , drop = FALSE]
          if (verbose) cat("Rows:", nrow(rows), " ")
          # If we have 0 rows for this bin in this fold, we need to debug.
          if (class(rows) != "data.frame" || nrow(rows) == 0) {
            #browser()
            NULL
          } else {
            rows
          }
        })

        # Remove elements that are NULL or 0 rows.
        for (element_i in length(combine_rows)) {
          item = combine_rows[[element_i]]
          if (class(item) != "data.frame" || nrow(item) == 0) {
            combine_rows[[element_i]] = NULL
          }
        }

        tryCatch({
        bin_df = do.call(rbind, combine_rows)
        if (verbose) cat("\n")

        # Create a list with one element ($val_preds df) per fold.
        bin_list = lapply(1:V, function(fold_i) {
          # Return with an enclosing list.
          list(bin_df[bin_df$fold == fold_i, ])
        })

        # Rename the element to be $val_preds
        for (fold in 1:V) {
          names(bin_list[[fold]]) = c("val_preds")
        }

        # bin_df can be NULL if the variable is skipped due to errors,
        # e.g. lack of variation.

        if (!is.null(bin_df) && nrow(bin_df) > 0L) {
          pooled_bin = estimate_pooled_results(bin_list, verbose = verbose)
          # Now we have $thetas and $influence_curves

          # Save the vector of estimates into the appropriate spot.
          # $thetas has the treatment-specific means
          # $influence_curves can be used to calculate the SE's, but shouldn't those
          # already be calculated by estimate_pooled_results()

          # Loop over fold results and insert the thetas into appropriate df.
          for (fold in 1:length(bin_list)) {
            bin_df = fold_results[[fold]]$bin_df
            row = bin_df$level == bin & bin_df$cv_fold == fold
            fold_results[[fold]]$bin_df[row, "test_theta_tmle"] = pooled_bin$thetas[fold]
            fold_results[[fold]]$bin_df[row, "test_var_tmle"] = var(pooled_bin$influence_curves[[fold]])
          }

        } else {
          cat("Skipping bin", bin, "- no rows are available.\n")
          # We have no val_preds for this bin, so skip pooled result estimation.

          # Temporary simplification for debugging purposes.
          #pooled_bin = list(thetas = 1:V)
          pooled_bin = list(thetas = rep(NA, V))
        }

        }, error = function(error) {
          pooled_bin = list(thetas = rep(NA, V))
        })


      }

        if (verbose) {
          cat("\n")
        }
      }

      # Combine results for each fold into a single dataframe.
      results_by_fold_and_level = do.call(rbind, lapply(fold_results, `[[`, "bin_df"))

      # Aggregate into a results_by_level dataframe.
      results_by_level = results_by_level(results_by_fold_and_level,
                                          verbose = verbose)

      ###########
      # Create list to save results for this variable.
      var_results = list(
        EY1V = NULL,
        EY0V = NULL,
        thetaV = NULL,
        varICV = NULL,
        labV = NULL,
        nV = NULL,
        fold_results = fold_results,
        type = "numeric",
        results_by_fold_and_level = results_by_fold_and_level,
        results_by_level = results_by_level,
        name = nameA
      )

      # TODO: compile results into the new estimate.

      pooled_min = estimate_pooled_results(lapply(fold_results, function(x) x$level_min),
                                           verbose = verbose)
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
            # Variance for the risk difference (maximal contrast parameter).
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
      # save(vim_numeric, file = "varimpact.RData")
      # TEMP remove this:
      stop(paste("We have", numerics$num_numeric, "continuous variables but only",
                 length(vim_numeric), "results."))
    }


    # Dataframe to hold all of the variable-by-fold-by-level results.
    #results_by_fold_and_level_obj = do.call(rbind, lapply(vim_numeric, `[[`, "results_by_fold_and_level"))

    # Dataframe to hold all of the variable-by-fold-by-level results.
    compile_results_by_fold_and_level = lapply(vim_numeric, function(result) {
      # Only extract results_by_fold_and_level if it's not NULL
      if ("results_by_fold_and_level" %in% names(result) &&
          !is.null(result$results_by_fold_and_level)# &&
          #!is.na(result$results_by_fold_and_level)
      ) {
        result$results_by_fold_and_level
      } else {
        NULL
      }
    })

    results_by_fold_and_level_obj = NULL
    # TODO: this shouldn't fail as easily - need to remove problematic results so that
    # the remainder can be used.
    tryCatch( {
      results_by_fold_and_level_obj = do.call(rbind, compile_results_by_fold_and_level)
    }, error = function(error) {
      cat("Errored while compiling results by fold and level.\n")
    })

    #results_by_level_obj = do.call(rbind, lapply(vim_numeric, `[[`, "results_by_level"))

    #results_by_level = do.call(rbind, lapply(vim_factor, `[[`, "results_by_level"))
    compile_results_by_level = lapply(vim_numeric, function(result) {
      # Only extract results_by_fold_and_level if it's not NULL
      if ("results_by_level" %in% names(result) &&
          !is.null(result$results_by_level) &&
          # TODO: figure out why expressions are going into this element.
          !is.expression(result$results_by_level)
          # !is.na(result$results_by_level)
      ) {
        result$results_by_level
      } else {
        NULL
      }
    })

    # TODO: this shouldn't fail as easily - need to remove problematic results so that
    # the remainder can be used.
    results_by_level_obj = NULL
    tryCatch({
      results_by_level_obj = do.call(rbind, compile_results_by_level)
    }, error = function(e) {
      # TODO: add browser?
      # TODO: figure out why this happens - presumably due to covariate that failed.
      # Error message:
      # Error in rep(xi, length.out = nvar) :
      #  attempt to replicate an object of type 'closure'
      cat("Errored while compiling results by level.\n")
    })

    colnames_numeric = colnames(numerics$data.cont.dist)
  } else {
    colnames_numeric = NULL
    vim_numeric = NULL
    results_by_fold_and_level_obj = NULL
    results_by_level_obj = NULL
    cat("No numeric variables for variable importance estimation.\n")
  }

  if (verbose) cat("Completed numeric variable importance estimation.\n")

  # Compile and return results.
  (results = list(
    vim_numeric = vim_numeric,
    results_by_fold_and_level = results_by_fold_and_level_obj,
    results_by_level = results_by_level_obj,
    colnames_numeric = colnames_numeric
    #, data.numW = numerics$data.numW
  ))
}
