test_that("Test whether ple_train works (ctns)", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  
  learners <- c("ranger", "glmnet", "linear", "bart")
  metas <- c("T-learner", "X-learner", "S-learner")
  res_dat <- NULL
  
  for (learn in learners) {
    for (meta in metas) {
      message(paste("Learner", learn))
      message(paste("Meta-learner", meta))
      res <- ple_train(Y, A, X, ple=learn, meta=meta)
      class_ind <- ifelse(class(res)=="ple_train", 1, 0)
      mu_ind <- ifelse(!is.null(res$mu_train), 1, 0)
      hold <- data.frame(learn=learn, meta=meta, class_ind=class_ind, mu_ind=mu_ind)
      res_dat <- rbind(res_dat, hold)
    }
  }
  res_ctns <- res_dat
  eql_ctns <- (mean(res_ctns$class_ind) + mean(res_ctns$mu_ind))/2
  expect_equal(eql_ctns, 1L)
})
test_that("Test whether ple_train works (binomial)", {
  
  skip_on_cran()
  ## Binomial ##
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  
  learners <- c("ranger", "glmnet", "linear", "bart")
  metas <- c("T-learner", "X-learner", "S-learner")
  res_dat <- NULL
  
  for (learn in learners) {
    for (meta in metas) {
      message(paste("Learner", learn))
      message(paste("Meta-learner", meta))
      res <- ple_train(Y, A, X, ple=learn, meta=meta, family="binomial")
      class_ind <- ifelse(class(res)=="ple_train", 1, 0)
      mu_ind <- ifelse(!is.null(res$mu_train), 1, 0)
      hold <- data.frame(learn=learn, meta=meta, class_ind=class_ind, mu_ind=mu_ind)
      res_dat <- rbind(res_dat, hold)
    }
  }
  res_bin <- res_dat
  eql_bin <- (mean(res_bin$class_ind) + mean(res_bin$mu_ind))/2
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
  
  learners <- c("ranger", "glmnet", "linear")
  metas <- c("T-learner", "S-learner")
  res_dat <- NULL
  
  for (learn in learners) {
    for (meta in metas) {
      message(paste("Learner", learn))
      message(paste("Meta-learner", meta))
      res <- ple_train(Y, A, X, ple=learn, meta=meta, family="survival")
      class_ind <- ifelse(class(res)=="ple_train", 1, 0)
      mu_ind <- ifelse(!is.null(res$mu_train), 1, 0)
      hold <- data.frame(learn=learn, meta=meta, class_ind=class_ind, mu_ind=mu_ind)
      res_dat <- rbind(res_dat, hold)
    }
  }
  res_surv <- res_dat
  eql_surv <- (mean(res_surv$class_ind) + mean(res_surv$mu_ind))/2
  expect_equal(eql_surv, 1L)
})