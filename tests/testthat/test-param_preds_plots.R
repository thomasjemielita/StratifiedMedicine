test_that("Test whether individual param models run alone and within PRISM; check forest plots", {
  
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
  param_checker = function(Y, A, X, ple=NULL, ple.hyper=NULL, filter="None", param=NULL,
                         submod=NULL, submod.hyper=NULL, family="gaussian",
                         verbose=FALSE){
    p0 <- PRISM(Y=Y, A=A, X=X, ple=ple, ple.hyper=ple.hyper, filter=filter, submod = submod,
                submod.hyper = submod.hyper, family=family, param=param,
                verbose = verbose)
    p1 <- do.call(param, list(Y, A, X, mu_hat=p0$mu_train, Subgrps=p0$out.train$Subgrps,
                              alpha_s = p0$alpha_s, alpha_ovrl = p0$alpha_ovrl,
                                     Xtest=X, family=family) )
    plt = plot(p0, type="forest")
    plt.test = ( "ggplot" %in% class(plt)   )
    return( data.frame(param=param, succeed=1, plt=plt.test))
  }
  
  # Loop through ple options #
  param.vec = c("param_lm", "param_ple", "param_dr")
  eql.dat = NULL
  for (param in param.vec){
    eqlz <- param_checker(Y=Y, A=A, X=Xstar,param=param)
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_ctns_default <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1] &
                                             sum(plt)==dim(eql.dat)[1],1,0) )
  
  ## Prognostic Test ##
  param.vec = c("param_lm", "param_ple")
  eql.dat = NULL
  for (param in param.vec){
    eqlz <- param_checker(Y=Y, A=NULL, X=Xstar,param=param)
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_ctns_prog <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1] &
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
  
  ## cox vs rmst ##
  param.vec = c("param_rmst", "param_cox")
  eql.dat = NULL
  for (param in param.vec){
    eqlz <- param_checker(Y=Y, A=A, X=Xstar,param=param, ple="None")
    eql.dat <- rbind(eql.dat, eqlz)
  }
  eql_surv_default <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1] &
                                          sum(plt)==dim(eql.dat)[1],1,0) )
  
  ## Prognostic Test ##
  eqlz <- param_checker(Y=Y, A=NULL, X=Xstar, ple="None", param="param_rmst")
  eql.dat <- eqlz
  eql_surv_prog <- with(eql.dat, ifelse(sum(succeed)==dim(eql.dat)[1] &
                                             sum(plt)==dim(eql.dat)[1],1,0) )
  
  
  ### Output Test Results ###
  expect_equal(eql_ctns_default, 1L)
  expect_equal(eql_ctns_prog, 1L)
  expect_equal(eql_surv_prog, 1L)
  expect_equal(eql_surv_prog, 1L)
})