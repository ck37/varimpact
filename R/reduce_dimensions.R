#' Reduce variables in a dataframe to a target number of covariates.
#'
#' Currently uses HOPACH hierarchical clustering but could be generalized.
#'
#' @param data Dataframe
#' @param newX Optional second dataframe to receive the same reduction.
#' @param max_variables Maximum we want to allow, after which dimension reduction
#'   will take place. Cannot be more than 15 due to HOPACH limitations. Set to NULL
#'   to disable any dimension reduction.
#' @param verbose If true will output additional information during execution.
#'
#' @importFrom hopach hopach distancematrix
#'
#' @export
reduce_dimensions = function(data, newX = NULL, max_variables, verbose = F) {

  num_columns = ncol(data)
  # Set this by default, then override it if we do reduce dimensions.
  variables = colnames(data)

  # Skip if number covariates is within the target maximum or if the maximum is null.
  if (num_columns <= max_variables || is.null(max_variables)) {
    Wtsht = data
    Wvsht = newX
  } else {
    if (verbose) cat("Reducing dimensions via clustering.\n")

    #mydist = as.matrix(hopach::distancematrix(t(Wt), d = "cosangle", na.rm = T))

    # Compute pairwise distances between each variable in the dataframe.
    # We transpose Wt because we want to cluster columns rather than rows.
    mydist = try(hopach::distancematrix(t(data), d = "cosangle", na.rm = T),
                 silent = !verbose)
    if (class(mydist) == "try-error") {
      cat("Error in HOPACH clustering: failed to calculate distance matrix.\n")
    }

    # Attempt #1.
    # We transpose Wt to cluster the columns rather than rows.
    # K = number of variables to choose.
    # kmax = maximum number of children at each node in the tree.
    # khigh = max # of children at each node when computing mss, usually the same.
    hopach.1 = try(hopach::hopach(t(data), dmat = mydist, mss = "mean", verbose = T,
                                  K = max_variables, kmax = 3, khigh = 3),
                   silent = !verbose)
    if (class(hopach.1) == "try-error") {
      if (verbose) {
        cat("Hopach attempt 1 fail.\n")
        print(hopach.1)
      }

      # Attempt #2.
      # We transpose Wt to cluster the columns rather than rows.
      hopach.1 <- try(hopach::hopach(t(data), dmat = mydist, mss = "med", verbose = T,
                                     K = max_variables, kmax = 3, khigh = 3),
                      silent = !verbose)
    }
    if (class(hopach.1) == "try-error") {
      if (verbose) {
        cat("Attempt 2 fail.")# Reverting to original W dataframe.\n")
        print(hopach.1)
      }

      # Attempt #3. Last try!
      # We transpose Wt to cluster the columns rather than rows.
      hopach.1 <- try(hopach::hopach(t(data), dmat = mydist, mss = "med", verbose = F,
                                     K = max_variables, kmax = 3, khigh = 3, newmed="nn"),
                      silent = !verbose)
    }
    if (class(hopach.1) == "try-error") {
      if (verbose) {
        cat("Attempt 3 fail. Reverting to original W dataframe.\n")
        # Now try to debug this.
        # stop()
      }
      #warning("Dimensionality reduction failed. i=", i, "V=", kk, "A=", nameA)
      Wtsht = data
      Wvsht = newX
    } else {
      # TODO: annotate what is going on here with the HOPACH result object.
      nlvls = nchar(max(hopach.1$final$labels))
      no = trunc(mean(log10(hopach.1$final$labels)))

      # Find highest level of tree where minimum number of covariates is >= adjust_cutoff.
      lvl = 1:nlvls
      ncv = NULL
      for (ii in lvl) {
        ncv = c(ncv, length(unique(trunc(hopach.1$final$labels/10^(no - (ii - 1))))))
      }
      ncv = unique(ncv)
      # Suppress possible "no non-missing arguments to min; returning Inf"
      # warning from min().
      # TODO: investigate more and confirm that this is ok.
      suppressWarnings({
        lev = min(min(nlvls, dim(data)[2]), min(lvl[ncv >= max_variables]))
      })
      two.clust <- unique(trunc(hopach.1$final$labels/(10^(no - (lev - 1)))))
      md <- hopach.1$final$medoids
      mm = md[, 1] %in% two.clust
      incc = md[mm, 2]
      Wtsht = data[, incc]
      Wvsht = newX[, incc]
      variables = colnames(data)[incc]
    }
    if (verbose) cat(" Updated columns:", ncol(Wtsht), "\n")
  }

  results = list(data = Wtsht, newX = Wvsht, variables = variables)
  results
}
