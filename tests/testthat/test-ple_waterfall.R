test_that("Test whether PLE waterfall plots are work", {
  
  skip_on_cran()
  
  ## Continuous ##
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  
  # Run Default: filter_glmnet, ple_ranger, lmtree, param_ple #
  res0 = PRISM(Y=Y, A=A, X=X)
  plot(res0, type="PLE:waterfall")
  
  res0 = ple_train(Y=Y, A=A, X=X)
  plot_ple(res0, type="waterfall")
  
  ## Binary ##
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  
  res0 = PRISM(Y=Y, A=A, X=X)
  plot(res0, type="PLE:waterfall")
  
  res0 = ple_train(Y=Y, A=A, X=X, family="binomial")
  plot_ple(res0, type="waterfall")
  
  
  
  library(survival)
  library(ggplot2)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom( n = dim(X)[1], size=1, prob=0.5  )
  
  # Linear #
  res1 = PRISM(Y=Y, A=A, X=X, ple = "linear")
  plot(res1, type="PLE:waterfall")
  
  res1 = ple_train(Y=Y, A=A, X=X, ple = "linear")
  plot_ple(res1, type="waterfall")
  
  # Glmnet #
  res2 = PRISM(Y=Y, A=A, X=X, ple = "glmnet", meta = "X-learner")
  plot(res2, type="PLE:waterfall")
  
  res2 = ple_train(Y=Y, A=A, X=X, ple = "glmnet")
  plot_ple(res2, type="waterfall")
  
})