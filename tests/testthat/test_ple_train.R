## Issues with binomial / ple_linear? ##
test_that("Test whether ple_train works (ctns)", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  
  # Should include bart? #
  learners <- c("ranger", "linear")
  metas <- c("T-learner", "X-learner", "S-learner", "DR-learner")
  res_dat <- NULL
  
  for (learn in learners) {
    for (meta in metas) {
      message(paste("Learner", learn))
      message(paste("Meta-learner", meta))
      res <- ple_train(Y, A, X, ple=learn, meta=meta)
      class_ind <- ifelse(class(res)=="ple_train", 1, 0)
      mu_ind <- ifelse(!is.null(res$mu_train), 1, 0)
      diff_ind <- ifelse(sum(is.na(res$mu_train$diff_1_0))==0, 1, 0)
      hold <- data.frame(learn=learn, meta=meta, 
                         class_ind=class_ind, mu_ind=mu_ind,
                         diff_ind = diff_ind)
      res_dat <- rbind(res_dat, hold)
    }
  }
  res_ctns <- res_dat
  eql_ctns <- (mean(res_ctns$class_ind) + mean(res_ctns$mu_ind) + mean(res_ctns$diff_ind))/3
  expect_equal(eql_ctns, 1L)
})
test_that("Test whether ple_train works (binomial)", {
  
  skip_on_cran()
  ## Binomial ##
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  
  # Should include bart? #
  learners <- c("ranger", "linear")
  metas <- c("T-learner", "X-learner", "S-learner", "DR-learner")
  res_dat <- NULL
  
  for (learn in learners) {
    for (meta in metas) {
      message(paste("Learner", learn))
      message(paste("Meta-learner", meta))
      res <- ple_train(Y, A, X, ple=learn, meta=meta, family="binomial")
      print(summary(res$mu_train))
      class_ind <- ifelse(class(res)=="ple_train", 1, 0)
      mu_ind <- ifelse(!is.null(res$mu_train), 1, 0)
      diff_ind <- ifelse(sum(is.na(res$mu_train$diff_1_0))==0, 1, 0)
      hold <- data.frame(learn=learn, meta=meta, class_ind=class_ind, 
                         mu_ind=mu_ind, diff_ind = diff_ind)
      res_dat <- rbind(res_dat, hold)
    }
  }
  res_bin <- res_dat
  eql_bin <- (mean(res_bin$class_ind) + mean(res_bin$mu_ind) + mean(res_bin$diff_ind))/3
  expect_equal(eql_bin, 1L)
})
test_that("Test whether ple_train works (survival)", {
  
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
  
  learners <- c("ranger", "linear")
  metas <- c("T-learner", "S-learner", "DR-learner")
  grid.dat <- expand.grid(learners, metas)
  tests <- !(grid.dat$Var1=="linear" & grid.dat$Var2=="DR-learner")
  grid.dat <- grid.dat[tests,]
  res_dat <- NULL
  
  for (ii in 1:dim(grid.dat)[1]) {
    learn <- grid.dat$Var1[ii]
    meta <- grid.dat$Var2[ii]
    message(paste("Learner", learn))
    message(paste("Meta-learner", meta))
    res <- ple_train(Y, A, X, ple=learn, meta=meta, family="survival")
    print(summary(res$mu_train))
    class_ind <- ifelse(class(res)=="ple_train", 1, 0)
    mu_ind <- ifelse(!is.null(res$mu_train), 1, 0)
    diff_ind <- ifelse(sum(is.na(res$mu_train$diff_1_0))==0, 1, 0)
    hold <- data.frame(learn=learn, meta=meta, class_ind=class_ind, 
                       mu_ind=mu_ind, diff_ind=diff_ind)
    res_dat <- rbind(res_dat, hold)
  }
  res_surv <- res_dat
  eql_surv <- (mean(res_surv$class_ind) + mean(res_surv$mu_ind) + mean(res_surv$diff_ind))/3
  expect_equal(eql_surv, 1L)
})