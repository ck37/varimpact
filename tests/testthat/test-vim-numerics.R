## Test 

test_that("adjustment set functionality works appropriately", {
  
          library(varimpact)
          library(SuperLearner)
          library(devtools)
          library(testthat)
          
          # Data setup taken from Chris Kennedy Github page (https://github.com/ck37/varimpact)
          set.seed(1, "L'Ecuyer-CMRG")
          N <- 300
          num_normal <- 5
          X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
          Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
          # Add some missing data to X so we can test imputation.
          for (i in 1:10) X[sample(nrow(X), 1), sample(ncol(X), 1)] <- NA
          
          Q_lib = c("SL.mean", "SL.glmnet", "SL.ranger", "SL.rpartPrune")
          g_lib = c("SL.mean", "SL.glmnet")
          
          #vim = varimpact(Y = Y, data = X, Q.library = Q_lib, g.library = g_lib)
          
          vim = varimpact(Y = Y, data = X, Q.library = Q_lib, g.library = g_lib, adjustment_exclusions = list("V1" = c("V2","V3"), "V3" = c("V1"), "V5" = c("V1","V2","V3")))
          
          # V1
          expect_equal(vim$all_vims$V1$fold_results[[1]]$bin_results[[1]]$W_names, c("V4","V5","Imiss_V1","Imiss_V2","Imiss_V3","Imiss_V5"))
          
          # V2
          expect_equal(vim$all_vims$V2$fold_results[[1]]$bin_results[[1]]$W_names, c("V1","V3","V4","V5","Imiss_V1","Imiss_V2","Imiss_V3","Imiss_V5"))
          
          # V3
          expect_equal(vim$all_vims$V3$fold_results[[1]]$bin_results[[1]]$W_names, c("V2","V4","V5","Imiss_V1","Imiss_V2","Imiss_V3","Imiss_V5"))
          
          # V4
          expect_equal(vim$all_vims$V4$fold_results[[1]]$bin_results[[1]]$W_names, c("V1","V2","V3","V5","Imiss_V1","Imiss_V2","Imiss_V3","Imiss_V5"))
          
          # V5
          expect_equal(vim$all_vims$V5$fold_results[[1]]$bin_results[[1]]$W_names, c("V4","Imiss_V1","Imiss_V2","Imiss_V3","Imiss_V5"))
          
          })



