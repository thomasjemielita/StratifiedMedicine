test_that("Test whether submod_train summary + pooling works (ctns)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  family <- "gaussian"
  
  # Run various configurations #
  res0 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "no")
  res1 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "no", 
                       R=5, resample = "Bootstrap")
  res2 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff", 
                       R = 5, resample = "Bootstrap")
  res3 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff_boot", 
                       R_pool = 5)
  res4 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff_boot", 
                       R_pool = 5, R = 2, resample = "Bootstrap")
  
  # Check #
  checker <- function(res, pool, submod) {
    
    summ <- summary(res)
    
    if (submod!="boot") {
      if (pool=="no") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (observed)`), 0, 1)
      }
      if (pool=="yes") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (optimal trt)`), 0, 1)
      }
    }
    if (submod=="boot") {
      if (pool=="no") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (observed)`), 0, 1)
      }
      if (pool=="yes") {
        ind0 <- ifelse(is.null(summ$`Treatment Effect Estimates (optimal trt)`), 0, 1)
        ind1 <- ifelse(is.null(summ$`Treatment Effect Estimates (bootstrap)`), 0, 1)
        ind <- mean(c(ind0, ind1))
      }
    }
    hold <- data.frame(pool = pool, submod=submod, ind = ind)
    return(hold)
  }
  
  out <- rbind(
    checker(res0, pool = "no", submod="no_boot"),
    checker(res1, pool = "no", submod="boot"),
    checker(res2, pool = "yes", submod="boot"),
    checker(res3, pool = "yes", submod="no_boot"),
    checker(res4, pool = "yes", submod="boot")
  )
  final_check <- ifelse(mean(out$ind)==1, TRUE, FALSE)
  
  
  # Check #
  expect_equal(final_check, TRUE)

})

test_that("Test whether submod_train summary + pooling works (binomial)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  library(pROC)
  
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  family <- "binomial"
  
  
  # Run various configurations #
  res0 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "no")
  res1 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "no", 
                       R=5, resample = "Bootstrap")
  res2 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff", 
                       R = 5, resample = "Bootstrap")
  res3 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff_boot", 
                       R_pool = 5)
  res4 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff_boot", 
                       R_pool = 5, R = 2, resample = "Bootstrap")
  
  # Check #
  checker <- function(res, pool, submod) {
    
    summ <- summary(res)
    
    if (submod!="boot") {
      if (pool=="no") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (observed)`), 0, 1)
      }
      if (pool=="yes") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (optimal trt)`), 0, 1)
      }
    }
    if (submod=="boot") {
      if (pool=="no") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (observed)`), 0, 1)
      }
      if (pool=="yes") {
        ind0 <- ifelse(is.null(summ$`Treatment Effect Estimates (optimal trt)`), 0, 1)
        ind1 <- ifelse(is.null(summ$`Treatment Effect Estimates (bootstrap)`), 0, 1)
        ind <- mean(c(ind0, ind1))
      }
    }
    hold <- data.frame(pool = pool, submod=submod, ind = ind)
    return(hold)
  }
  
  out <- rbind(
    checker(res0, pool = "no", submod="no_boot"),
    checker(res1, pool = "no", submod="boot"),
    checker(res2, pool = "yes", submod="boot"),
    checker(res3, pool = "yes", submod="no_boot"),
    checker(res4, pool = "yes", submod="boot")
  )
  final_check <- ifelse(mean(out$ind)==1, TRUE, FALSE)
  
  
  # Check #
  expect_equal(final_check, TRUE)
  
})

test_that("Test whether submod_train summary + pooling works (survival)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  library(survival)
  require(TH.data); require(coin)
  library(pROC)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom(n = dim(X)[1], size=1, prob=0.5)
  
  family <- "survival"
  
  # Run various configurations #
  res0 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "no")
  res1 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "no", 
                       R=5, resample = "Bootstrap")
  res2 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff", 
                       R = 5, resample = "Bootstrap")
  res3 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff_boot", 
                       R_pool = 5)
  res4 <- submod_train(Y=Y, A=A, X=X, submod = "lmtree", pool = "trteff_boot", 
                       R_pool = 5, R = 2, resample = "Bootstrap")
  
  # Check #
  checker <- function(res, pool, submod) {
    
    summ <- summary(res)
    
    if (submod!="boot") {
      if (pool=="no") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (observed)`), 0, 1)
      }
      if (pool=="yes") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (optimal trt)`), 0, 1)
      }
    }
    if (submod=="boot") {
      if (pool=="no") {
        ind <- ifelse(is.null(summ$`Treatment Effect Estimates (observed)`), 0, 1)
      }
      if (pool=="yes") {
        ind0 <- ifelse(is.null(summ$`Treatment Effect Estimates (optimal trt)`), 0, 1)
        ind1 <- ifelse(is.null(summ$`Treatment Effect Estimates (bootstrap)`), 0, 1)
        ind <- mean(c(ind0, ind1))
      }
    }
    hold <- data.frame(pool = pool, submod=submod, ind = ind)
    return(hold)
  }
  
  out <- rbind(
    checker(res0, pool = "no", submod="no_boot"),
    checker(res1, pool = "no", submod="boot"),
    checker(res2, pool = "yes", submod="boot"),
    checker(res3, pool = "yes", submod="no_boot"),
    checker(res4, pool = "yes", submod="boot")
  )
  final_check <- ifelse(mean(out$ind)==1, TRUE, FALSE)
  
  
  # Check #
  expect_equal(final_check, TRUE)
 
  
})
