test_that("Test whether submod models match PRISM; check plot", {
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  # Filter #
  fmod = filter_glmnet(Y, A, X)
  Xstar = X[,colnames(X) %in% fmod$filter.vars]
  
  ## submod checker (does PRISM output match wrapper output; can we get ggplot?) ##
  submod_checker = function(Y, A, X, ple=NULL, filter="None", 
                            submod, submod.hyper=NULL, 
                            verbose=FALSE){
    p0 <- PRISM(Y=Y, A=A, X=X, ple=ple, filter=filter, submod = submod,
                submod.hyper = submod.hyper,
                verbose = verbose)
    p1 <- do.call(submod, append(list(Y, A, X, Xtest=X, mu_train=p0$mu_train),
                                 submod.hyper) )
    eql_1a = all.equal(p1, p0$submod.fit)
    plt_1 = plot(p0)
    eql_1b = ( "ggplot" %in% class(plt_1)   )
    eql.dat = data.frame(submod=submod, fit=eql_1a, plt = eql_1b)
    return(eql.dat)
  }

  # Loop through submod options #
  submod.vec = c("submod_lmtree", "submod_ctree", "submod_otr", "submod_rpart")
  eql.dat = NULL
  for (submod in submod.vec){
    eqlz <- submod_checker(Y=Y, A=A, X=Xstar, submod=submod)
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_ctns_default <- with(eql.dat, ifelse(sum(fit)==dim(eql.dat)[1] & 
                                   sum(plt)==dim(eql.dat)[1],1,0) )
  
  ## Prognostic Test ##
  eql_ctns_prog <- submod_checker(Y=Y, A=NULL, X=Xstar, submod="submod_ctree")
  eql_ctns_prog <- with(eql_ctns_prog, ifelse(fit==TRUE & plt==TRUE, 1, 0))
  
  ## PLE versions (ctree, rpart, otr) ##
  submod.vec = c("submod_ctree", "submod_otr", "submod_rpart")
  eql.dat = NULL
  for (submod in submod.vec){
    eqlz <- submod_checker(Y=Y, A=A, X=Xstar, submod=submod, 
                           submod.hyper = list(outcome_PLE=TRUE))
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_ctns_ple <- with(eql.dat, ifelse(sum(fit)==dim(eql.dat)[1] & 
                                             sum(plt)==dim(eql.dat)[1],1,0) )
  
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
  submod.vec = c("submod_ctree", "submod_weibull")
  eql.dat = NULL
  for (submod in submod.vec){
    eqlz <- submod_checker(Y=Y, A=A, X=Xstar, submod=submod, ple="None")
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_surv_default <- with(eql.dat, ifelse(sum(fit)==dim(eql.dat)[1] & 
                                         sum(plt)==dim(eql.dat)[1],1,0) )
  
  ## Prognostic Test ##
  eql_surv_prog <- submod_checker(Y=Y, A=NULL, X=Xstar, ple="None",
                                  submod="submod_ctree")
  eql_surv_prog <- with(eql_surv_prog, ifelse(fit==TRUE & plt==TRUE, 1, 0))
  
  
  ### Output Test Results ###
  expect_equal(eql_ctns_default, 1L)
  expect_equal(eql_ctns_prog, 1L)
  expect_equal(eql_ctns_ple, 1L)
  expect_equal(eql_surv_prog, 1L)
  expect_equal(eql_surv_prog, 1L)
})