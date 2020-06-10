test_that("Test whether param_est works (ctns)", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  
  # Fit Subgroup Model #
  res_s <- submod_train(Y, A, X, submod="lmtree")
  Subgrps <- res_s$Subgrps.train
  
  # PLE #
  res_p <- ple_train(Y, A, X)
  mu_hat <- res_p$mu_train
  
  params <- c("lm", "dr", "ple")
  res_dat <- NULL
  
  for (param in params) {
    message(paste("Parameter Estimation", param))
    res <- param_est(Y, A, X, param=param, mu_hat=mu_hat,
                     Subgrps=Subgrps)
    data_ind <- ifelse(is.data.frame(res), 1, 0)
    hold <- data.frame(param=param, data_ind=data_ind)
    res_dat <- rbind(res_dat, hold)
  }
  res_ctns <- res_dat
  eql_ctns <- (mean(res_ctns$data_ind))
  expect_equal(eql_ctns, 1L)
})
test_that("Test whether ple_train works (binomial)", {
  
  skip_on_cran()
  ## Binomial ##
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  # Fit Subgroup Model #
  res_s <- submod_train(Y, A, X, submod="glmtree")
  Subgrps <- res_s$Subgrps.train
  
  # PLE #
  res_p <- ple_train(Y, A, X)
  mu_hat <- res_p$mu_train
  
  params <- c("lm", "dr", "ple")
  res_dat <- NULL
  
  for (param in params) {
    message(paste("Parameter Estimation", param))
    res <- param_est(Y, A, X, param=param, mu_hat=mu_hat,
                     Subgrps=Subgrps)
    data_ind <- ifelse(is.data.frame(res), 1, 0)
    hold <- data.frame(param=param, data_ind=data_ind)
    res_dat <- rbind(res_dat, hold)
  }
  res_bin <- res_dat
  eql_bin <- (mean(res_bin$data_ind))
  expect_equal(eql_bin, 1L)
})
test_that("Test whether ple_train works (binomial)", {
  
  skip_on_cran()
  ### Survival Tests ###
  library(survival)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom(n = dim(X)[1], size=1, prob=0.5)
  
  # Fit Subgroup Model #
  res_s <- submod_train(Y, A, X, submod="mob_weib")
  Subgrps <- res_s$Subgrps.train
  
  # PLE #
  res_p <- ple_train(Y, A, X)
  mu_hat <- res_p$mu_train
  
  params <- c("cox", "rmst")
  res_dat <- NULL
  
  for (param in params) {
    message(paste("Parameter Estimation", param))
    res <- param_est(Y, A, X, param=param, mu_hat=mu_hat,
                     Subgrps=Subgrps)
    data_ind <- ifelse(is.data.frame(res), 1, 0)
    hold <- data.frame(param=param, data_ind=data_ind)
    res_dat <- rbind(res_dat, hold)
  }
  res_surv <- res_dat
  eql_surv <- (mean(res_surv$data_ind))
  expect_equal(eql_surv, 1L)
})