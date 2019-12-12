test_that("Test whether individual ple models run alone and within PRISM; check plots", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  # Filter #
  fmod = filter_glmnet(Y, A, X)
  Xstar = X[,colnames(X) %in% fmod$filter.vars]
  
  ## ple checker (does PRISM output match wrapper output; can we get ggplot?) ##
  ple_checker = function(Y, A, X, ple=NULL, ple.hyper=NULL, filter="None", 
                            submod=NULL, submod.hyper=NULL, family="gaussian",
                            verbose=FALSE){
    p0 <- PRISM(Y=Y, A=A, X=X, ple=ple, ple.hyper=ple.hyper, filter=filter, submod = submod,
                submod.hyper = submod.hyper, family=family,
                verbose = verbose)
    p1 <- do.call(ple, append(list(Y, A, X, Xtest=X, family=family),
                                 ple.hyper) )
    plt_1 = plot(p0, type="PLE:density")
    return( data.frame(ple=ple, succeed=1))
  }
  
  # Loop through ple options #
  ple.vec = c("ple_ranger", "ple_bart", "ple_glmnet")
  eql.dat = NULL
  for (ple in ple.vec){
    eqlz <- ple_checker(Y=Y, A=A, X=Xstar,ple=ple)
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_ctns_default <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1],1,0) )
  
  ## Prognostic Test ##
  ple.vec = c("ple_ranger", "ple_bart", "ple_glmnet")
  eql.dat = NULL
  for (ple in ple.vec){
    eqlz <- ple_checker(Y=Y, A=NULL, X=Xstar,ple=ple)
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_ctns_prog <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1],1,0) )
  
  ## Try a few more variations ##
  hold1 <- ple_checker(Y=Y, A=A, X=Xstar,ple="ple_ranger", ple.hyper=list(byTrt=FALSE))
  eql_ctns_hyper = ifelse(sum(hold1$succeed)==dim(hold1)[1],1,0)
  
  ### Survival Tests ###
  library(survival)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
  # Filter #
  fmod = filter_glmnet(Y, A, X, family="survival")
  Xstar = X[,colnames(X) %in% fmod$filter.vars]
  
  ## Standard (weibull, ctree) ##
  ple.vec = c("ple_ranger", "ple_glmnet")
  eql.dat = NULL
  for (ple in ple.vec){
    eqlz <- ple_checker(Y=Y, A=A, X=Xstar,ple=ple, family="survival")
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_surv_default <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1],1,0) )
  
  ## Prognostic Test ##
  for (ple in ple.vec){
    eqlz <- ple_checker(Y=Y, A=NULL, X=Xstar,ple=ple, family="survival")
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_surv_prog <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1],1,0) )
  
  
  ### Output Test Results ###
  expect_equal(eql_ctns_default, 1L)
  expect_equal(eql_ctns_prog, 1L)
  expect_equal(eql_ctns_hyper, 1L)
  expect_equal(eql_surv_prog, 1L)
  expect_equal(eql_surv_prog, 1L)
})