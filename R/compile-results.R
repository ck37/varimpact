compile_results =
  function(colnames_numeric,
           colnames_factor,
           vim_numeric,
           vim_factor,
           V,
           verbose = FALSE) {

  num_numeric = length(colnames_numeric)
  num_factor = length(colnames_factor)

  variable_types = c(rep("Ordered", num_numeric), rep("Factor", num_factor))
  variable_names = c(colnames_numeric, colnames_factor)

  all_vims = c(vim_numeric, vim_factor)
  names(all_vims) = variable_names

  element_length = sapply(all_vims, length)

  # May increase this to 8, hence >= operator in following lines.
  expected_length = 8

  num_short_vars = sum(element_length < expected_length)
  if (verbose && num_short_vars > 0) {
    # Identify which variables do not have expected length.
    cat(num_short_vars, "variables did not contain >=", expected_length,
        "result elements. Removing those variables from the results.\n")
    # Likely missing EY1V, EY0V, labV, and fold_results.
    print(element_length[element_length < expected_length])
    #     browser()
  }

  # Restrict to vim results that have at least 7 elements.
  vim_combined = all_vims[element_length >= expected_length]
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
      #cat("Lengths:", element_length, "\n")
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

      # Get labels for each of the training samples.
      extract_labels = function(x, total_folds) {
        # Repeat each fold id twice, one for low level and once for high.
        labels = rep(1:total_folds, 2)
        # Re-order in ascending order.
        # TODO: isn't there another rep() function that would sort automatically?
        oo = order(labels)
        labels = labels[oo]
        # Convert labels from vector elements to columns in a single row.
        out = as.vector(t(x))
        names(out) = paste0(rep(c("Low", "High"), total_folds), "_v", labels)
        out
      }

      # labV is result element 5.
      tst = lapply(results_no_na, function(x) x$labV)
      tst = lapply(tst, extract_labels, total_folds = V)
      labels = do.call(rbind, tst)

      # Each row is a variable and each column in a fold estimate.
      meanvarIC = apply(varIC, 1, mean)

      psi = apply(theta, 1, mean)
      SE = sqrt(meanvarIC / n)

      ci_lower = psi - 1.96 * SE
      ci_upper = psi + 1.96 * SE

      # Number of significant digits.
      signif_digits = 3

      # TODO: provide ci_lower and ci_upper as separate elements.
      CI95 = paste0("(", signif(ci_lower, signif_digits), " - ", signif(ci_upper, signif_digits), ")")

      # 1-sided p-value
      pvalue = 1 - pnorm(psi / SE)

      ##### FOR THETA (generalize to chi-square test?)
      # TT = (theta[,1] - theta[,2]) / sqrt(SEV[,1]^2 + SEV[,2]^2)
      # pval.comp=2*(1-pnorm(abs(TT))) FOR levels
      # (just make sure in same order)

      num_continuous = sum(variable_types == "Ordered")
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

        # For numeric variables, consistency means each fold finds the same directionality.
        cons = apply(dir, 1, length.uniq)
      }

      ##################
      # Factors
      num_factors = num_vars - num_continuous
      if (num_factors > 0) {
        lwr = NULL
        uwr = NULL
        for (i in 1:V) {
          # The 2 here is because we have the a_l and a_h labels, not because V = 2.
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
                            theta[sorted_rows, , drop = FALSE],
                            psi[sorted_rows],
                            CI95[sorted_rows],
                            res$adj,
                            labels[sorted_rows, , drop = FALSE],
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
                            labels,
                            consist)
      } else {
        outres = NULL
      }

      # TODO: this will give an error if we have no results.

      # Restrict to variables that aren't missing their p-value.
      outres = outres[!is.na(outres[, "rawp"]), , drop = F]
      #print(colnames(outres))
      #print(ncol(outres))

      #names(outres)[1:(1 + 2*V)] = c("VarType", paste0("Est_v", 1:V), "AvePsi", "CI95")
      names(outres)[1:(3 + V)] = c("VarType", paste0("Est_v", 1:V), "AvePsi", "CI95")
      #names(outres)[(9 + 2 * V)] = "Consistent"
      names(outres)[ncol(outres)] = "Consistent"

      #print(colnames(outres))
      #print(ncol(outres))

      ################
      # Get Consistency Measure and only significant
      # TODO: Make BH cut-off flexible in future versions (default at 0.05)
      outres.cons = outres[outres[, "BH"] < 0.05, , drop = F]
      outres.cons = outres.cons[outres.cons[, "Consistent"],
                                c("VarType", "AvePsi", "rawp", "BH", "CI95"), drop = F]
      colnames(outres.cons) = c("Type", "Estimate", "P-value", "Adj. P-value", "CI 95")


      # drops = c('VarType','description','Holm,')
      # outres.all=outres[,!(names(outres) %in% drops)]
      # We want to extract the per-fold psi estimates, per-fold levels, and consistency flag.
      outres.byV = outres[, c(2:(2 + V - 1), (7 + V):(7 + 3 * V)), drop = F]
      outres.all = outres[, c("VarType", "AvePsi", "CI95", "rawp", "BH", "Consistent"), drop = F]
      colnames(outres.all) = c("Type", "Estimate", "CI95", "P-value", "Adj. p-value", "Consistent")
    }
  }

  # Return results.
  results = list(results_consistent = outres.cons,
                 results_all = outres.all,
                 results_by_fold = outres.byV,
                 results_raw = outres,
                 all_vims = all_vims)


  return(results)
}
