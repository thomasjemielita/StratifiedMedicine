test_that("Test whether filter_train works (gaussian, binomial, survival)", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  
  ## Glmnet ##
  mod1 <- filter_train(Y, A, X, filter="glmnet")
  ind1 <- ifelse(class(mod1)=="filter_train", 1, 0)
  mod2 <- filter_train(Y, A, X, filter="glmnet", hyper = list(interaction=T))
  ind2 <- ifelse(class(mod2)=="filter_train", 1, 0)
  ## Ranger ##
  mod3 <- filter_train(Y, A, X, filter="ranger")
  ind3 <- ifelse(class(mod3)=="filter_train", 1, 0)
  sum_class <- sum(ind1, ind2, ind3)
  # Check if filter variables are in original covariate space #
  var1 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var2 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var3 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var_lengths <- var1+var2+var3
  eql_ctns <- ifelse(sum_class==3 & var_lengths==0,1,0)
  
  ## Binomial ##
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  
  ## Glmnet ##
  mod1 <- filter_train(Y, A, X, filter="glmnet", family="binomial")
  ind1 <- ifelse(class(mod1)=="filter_train", 1, 0)
  mod2 <- filter_train(Y, A, X, filter="glmnet", hyper = list(interaction=T),
                       family="binomial")
  ind2 <- ifelse(class(mod2)=="filter_train", 1, 0)
  ## Ranger ##
  mod3 <- filter_train(Y, A, X, filter="ranger", family="binomial")
  ind3 <- ifelse(class(mod3)=="filter_train", 1, 0)
  sum_class <- sum(ind1, ind2, ind3)
  # Check if filter variables are in original covariate space #
  var1 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var2 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var3 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var_lengths <- var1+var2+var3
  eql_bin <- ifelse(sum_class==3 & var_lengths==0,1,0)
  
  ### Survival Tests ###
  library(survival)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom( n = dim(X)[1], size=1, prob=0.5)
  
  ## Glmnet ##
  mod1 <- filter_train(Y, A, X, filter="glmnet")
  ind1 <- ifelse(class(mod1)=="filter_train", 1, 0)
  mod2 <- filter_train(Y, A, X, filter="glmnet", hyper = list(interaction=T))
  ind2 <- ifelse(class(mod2)=="filter_train", 1, 0)
  ## Ranger ##
  mod3 <- filter_train(Y, A, X, filter="ranger")
  ind3 <- ifelse(class(mod3)=="filter_train", 1, 0)
  sum_class <- sum(ind1, ind2, ind3)
  # Check if filter variables are in original covariate space #
  var1 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var2 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var3 <- length(setdiff(mod1$filter.vars, colnames(X)))
  var_lengths <- var1+var2+var3
  eql_surv <- ifelse(sum_class==3 & var_lengths==0,1,0)
  
  ### Output Test Results ###
  expect_equal(eql_ctns, 1L)
  expect_equal(eql_bin, 1L)
  expect_equal(eql_surv, 1L)
})