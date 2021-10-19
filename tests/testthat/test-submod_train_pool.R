test_that("Test whether submod_train pooling works (ctns)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  library(pROC)
  
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  family <- "gaussian"
  
  # Estimate Counterfactuals #
  ple_mod <- ple_train(Y=Y, A=A, X=X, family = family)
  mu_train <- ple_mod$mu_train
  
  delta_num <- 0.05
  delta <- "> 0.05"
  
  # Fit lmtree without submod #
  sub1 <- lmtree(Y~A | . , data = X, alpha = 0.10, 
                 maxdepth = 3, minsize = floor(dim(X)[1]*0.10))
  
  pred_s1 <- as.numeric(predict(sub1, type="node", newdata = X))
  
  param_s1 <- param_est(Y=Y, A=A, X=X, Subgrps = pred_s1, param="lm")
  param_s1$dopt <- ifelse(param_s1$est > delta_num , "dopt=1", "dopt=0")
  
  trt_assign0 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign0 <- dplyr::left_join(trt_assign0, param_s1, by="Subgrps")
  
  # Pooling (trteff) #
  mod0 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="lmtree", pool = "trteff", 
                       hyper = list(maxdepth=3), delta = delta)
  
  eql_trteff <- all.equal(trt_assign0[,c("Subgrps", "dopt")], 
                          mod0$trt_assign)
  
  # Set up data #
  ind_ple <- ifelse(mu_train$diff_1_0 > delta_num, 1, 0)
  w_ple <- abs(mu_train$diff_1_0 - delta_num)
  Subgrps.mat <- data.frame(Subgrps=as.factor(pred_s1))
  Subgrps.mat <- data.frame(model.matrix(~.-1, data=Subgrps.mat))
  otr_dat <- data.frame(ind_ple, Subgrps.mat)
  
  # logistic (only subgroups) #
  mod <- suppressWarnings(glm(ind_ple ~ . -1, 
                                data = otr_dat,
                                family = "binomial", weights = w_ple))
  prob_dopt1 <- as.numeric(predict(mod, type="response"))
  prob_dopt1 <- aggregate(prob_dopt1 ~ pred_s1, FUN="mean")
  colnames(prob_dopt1) <- c("Subgrps", "prob_dopt")
  prob_dopt1$Subgrps <- as.character(prob_dopt1$Subgrps)
  
  trt_assign1 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign1 <- dplyr::left_join(trt_assign1, prob_dopt1, by="Subgrps")
  
  rocobj <- suppressMessages(pROC::roc(ind_ple, trt_assign1$prob_dopt))
  yindex1 <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  trt_assign1$dopt <- with(trt_assign1, ifelse(prob_dopt <= yindex1, "dopt=0",
                                               "dopt=1"))

  
  # ranger: include X again #
  mod_rf <- ranger::ranger(ind_ple ~ ., 
                          data=data.frame(otr_dat, X),
                          case.weights = w_ple)
  prob_dopt2 <- mod_rf$predictions
  prob_dopt2 <- aggregate(prob_dopt2 ~ pred_s1, FUN="mean")
  colnames(prob_dopt2) <- c("Subgrps", "prob_dopt")
  prob_dopt2$Subgrps <- as.character(prob_dopt2$Subgrps)
  
  trt_assign2 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign2 <- dplyr::left_join(trt_assign2, prob_dopt2, by=c("Subgrps"))
  
  rocobj <- suppressMessages(pROC::roc(ind_ple, trt_assign2$prob_dopt))
  yindex2 <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  trt_assign2$dopt <- with(trt_assign2, ifelse(prob_dopt <= yindex2, "dopt=0",
                                               "dopt=1"))
  
  # Run SUBMOD #
  mod1 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="lmtree", pool = "otr:logistic", 
                       hyper = list(maxdepth=3), delta = delta,
                       otr_prob_thres = "youden")
  
  eql_otrlog <- all.equal(trt_assign1[,c("Subgrps", "dopt")], 
                          mod1$trt_assign)
  
  mod2 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="lmtree", pool = "otr:rf", 
                       hyper = list(maxdepth=3), delta = delta,
                       otr_prob_thres = "youden")
  
  eql_otrrf <- all.equal(trt_assign2[,c("Subgrps", "dopt")], 
                          mod2$trt_assign)
  
  # Check if the youden index's are similar #
  diff_yindex1 <- abs(yindex1 - unique(mod1$trt_eff_pool$thres_youden))
  diff_yindex2 <- abs(yindex2 - unique(mod2$trt_eff_pool$thres_youden))
  
  bound <- 0.05
  diff_yindex1 <- ifelse(diff_yindex1 < bound, TRUE, FALSE)
  diff_yindex2 <- ifelse(diff_yindex2 < bound, TRUE, FALSE)
  
  # Check #
  expect_equal(eql_trteff, TRUE)
  expect_equal(eql_otrlog, TRUE)
  expect_equal(eql_otrrf, TRUE)
  expect_equal(diff_yindex1, TRUE)
  expect_equal(diff_yindex2, TRUE)
  
})

test_that("Test whether submod_train pooling works (binomial)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  library(pROC)
  
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  family <- "binomial"
  

  # Estimate Counterfactuals #
  ple_mod <- ple_train(Y=Y, A=A, X=X, family = family)
  mu_train <- ple_mod$mu_train
  
  delta_num <- 0.05
  delta <- "> 0.05"
  
  # Fit mob without submod #
  sub1 <- glmtree(Y~A | . , data = X, alpha = 0.10, 
                  family= binomial(link="identity"), 
                  maxdepth = 3, minsize = floor(dim(X)[1]*0.10))
  
  pred_s1 <- as.numeric(predict(sub1, type="node", newdata = X))
  
  param_s1 <- param_est(Y=Y, A=A, X=X, Subgrps = pred_s1, param="lm")
  param_s1$dopt <- ifelse(param_s1$est > delta_num , "dopt=1", "dopt=0")
  
  trt_assign0 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign0 <- dplyr::left_join(trt_assign0, param_s1, by="Subgrps")
  
  # Pooling (trteff) #
  mod0 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="glmtree", pool = "trteff", 
                       hyper = list(maxdepth=3), delta = delta)
  
  eql_trteff <- all.equal(trt_assign0[,c("Subgrps", "dopt")], 
                          mod0$trt_assign)
  
  # Set up data #
  ind_ple <- ifelse(mu_train$diff_1_0 > delta_num, 1, 0)
  w_ple <- abs(mu_train$diff_1_0 - delta_num)
  Subgrps.mat <- data.frame(Subgrps=as.factor(pred_s1))
  Subgrps.mat <- data.frame(model.matrix(~.-1, data=Subgrps.mat))
  otr_dat <- data.frame(ind_ple, Subgrps.mat)
  
  # logistic (only subgroups) #
  mod <- suppressWarnings(glm(ind_ple ~ . -1, 
                              data = otr_dat,
                              family = "binomial", weights = w_ple))
  prob_dopt1 <- as.numeric(predict(mod, type="response"))
  prob_dopt1 <- aggregate(prob_dopt1 ~ pred_s1, FUN="mean")
  colnames(prob_dopt1) <- c("Subgrps", "prob_dopt")
  prob_dopt1$Subgrps <- as.character(prob_dopt1$Subgrps)
  
  trt_assign1 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign1 <- dplyr::left_join(trt_assign1, prob_dopt1, by="Subgrps")
  
  rocobj <- suppressMessages(pROC::roc(ind_ple, trt_assign1$prob_dopt))
  yindex1 <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  trt_assign1$dopt <- with(trt_assign1, ifelse(prob_dopt <= yindex1, "dopt=0",
                                               "dopt=1"))
  
  
  # ranger: include X again #
  mod_rf <- ranger::ranger(ind_ple ~ ., 
                           data=data.frame(otr_dat, X),
                           case.weights = w_ple)
  prob_dopt2 <- mod_rf$predictions
  prob_dopt2 <- aggregate(prob_dopt2 ~ pred_s1, FUN="mean")
  colnames(prob_dopt2) <- c("Subgrps", "prob_dopt")
  prob_dopt2$Subgrps <- as.character(prob_dopt2$Subgrps)
  
  trt_assign2 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign2 <- dplyr::left_join(trt_assign2, prob_dopt2, by=c("Subgrps"))
  
  rocobj <- suppressMessages(pROC::roc(ind_ple, trt_assign2$prob_dopt))
  yindex2 <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  trt_assign2$dopt <- with(trt_assign2, ifelse(prob_dopt <= yindex2, "dopt=0",
                                               "dopt=1"))
  
  # Run SUBMOD #
  mod1 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="glmtree", pool = "otr:logistic", 
                       hyper = list(maxdepth=3), delta = delta,
                       otr_prob_thres = "youden")
  
  eql_otrlog <- all.equal(trt_assign1[,c("Subgrps", "dopt")], 
                          mod1$trt_assign)
  
  mod2 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="glmtree", pool = "otr:rf", 
                       hyper = list(maxdepth=3), delta = delta,
                       otr_prob_thres = "youden")
  
  eql_otrrf <- all.equal(trt_assign2[,c("Subgrps", "dopt")], 
                         mod2$trt_assign)
  
  # Check if the youden index's are similar #
  diff_yindex1 <- abs(yindex1 - unique(mod1$trt_eff_pool$thres_youden))
  diff_yindex2 <- abs(yindex2 - unique(mod2$trt_eff_pool$thres_youden))
  
  bound <- 0.05
  diff_yindex1 <- ifelse(diff_yindex1 < bound, TRUE, FALSE)
  diff_yindex2 <- ifelse(diff_yindex2 < bound, TRUE, FALSE)
  
  # Check #
  expect_equal(eql_trteff, TRUE)
  expect_equal(eql_otrlog, TRUE)
  expect_equal(eql_otrrf, TRUE)
  expect_equal(diff_yindex1, TRUE)
  expect_equal(diff_yindex2, TRUE)
  
})

test_that("Test whether submod_train pooling works (survival)", {
  
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
  
  # Estimate Counterfactuals #
  ple_mod <- ple_train(Y=Y, A=A, X=X, family = family)
  mu_train <- ple_mod$mu_train
  
  # Fit mob without submod #
  # MOB-WEIB #
  wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...) {
    survreg(y ~ 0 + x, weights = weights, dist = "weibull", ...)
  }
  logLik.survreg <- function(object, ...) {
    structure(object$loglik[2], df = sum(object$df), class = "logLik")
  }
  sub1 <- partykit::mob(Y ~ A | ., data = X,
                        fit = wbreg, 
                        control = partykit::mob_control(alpha=0.10, 
                                                        minsize=floor(dim(X)[1]*0.10),
                                                        maxdepth=3))
  
  delta_num <- 0
  delta <- "> 0"
  
  pred_s1 <- as.numeric(predict(sub1, type="node", newdata = X))
  
  param_s1 <- param_est(Y=Y, A=A, X=X, Subgrps = pred_s1, param="aft")
  param_s1$dopt <- ifelse(param_s1$est > delta_num , "dopt=1", "dopt=0")
  
  trt_assign0 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign0 <- dplyr::left_join(trt_assign0, param_s1, by="Subgrps")
  
  # Pooling (trteff) #
  mod0 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train, param = "aft",
                       submod="mob_weib", pool = "trteff", 
                       hyper = list(maxdepth=3), delta = delta)
  
  eql_trteff <- all.equal(trt_assign0[,c("Subgrps", "dopt")], 
                          mod0$trt_assign)
  
  # Set up data #
  delta_num <- 0.11
  delta <- "> 0.11"
  
  ind_ple <- ifelse(mu_train$diff_1_0 > delta_num, 1, 0)
  w_ple <- abs(mu_train$diff_1_0 - delta_num)
  Subgrps.mat <- data.frame(Subgrps=as.factor(pred_s1))
  Subgrps.mat <- data.frame(model.matrix(~.-1, data=Subgrps.mat))
  otr_dat <- data.frame(ind_ple, Subgrps.mat)
  
  # logistic (only subgroups) #
  mod <- suppressWarnings(glm(ind_ple ~ . -1, 
                              data = otr_dat,
                              family = "binomial", weights = w_ple))
  prob_dopt1 <- as.numeric(predict(mod, type="response"))
  prob_dopt1 <- aggregate(prob_dopt1 ~ pred_s1, FUN="mean")
  colnames(prob_dopt1) <- c("Subgrps", "prob_dopt")
  prob_dopt1$Subgrps <- as.character(prob_dopt1$Subgrps)
  
  trt_assign1 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign1 <- dplyr::left_join(trt_assign1, prob_dopt1, by="Subgrps")
  
  rocobj <- suppressMessages(pROC::roc(ind_ple, trt_assign1$prob_dopt))
  yindex1 <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  trt_assign1$dopt <- with(trt_assign1, ifelse(prob_dopt <= yindex1, "dopt=0",
                                               "dopt=1"))
  
  
  # ranger: include X again #
  mod_rf <- ranger::ranger(ind_ple ~ ., 
                           data=data.frame(otr_dat, X),
                           case.weights = w_ple)
  prob_dopt2 <- mod_rf$predictions
  prob_dopt2 <- aggregate(prob_dopt2 ~ pred_s1, FUN="mean")
  colnames(prob_dopt2) <- c("Subgrps", "prob_dopt")
  prob_dopt2$Subgrps <- as.character(prob_dopt2$Subgrps)
  
  trt_assign2 <- data.frame(Subgrps = as.character(pred_s1))
  trt_assign2 <- dplyr::left_join(trt_assign2, prob_dopt2, by=c("Subgrps"))
  
  rocobj <- suppressMessages(pROC::roc(ind_ple, trt_assign2$prob_dopt))
  yindex2 <- suppressWarnings(pROC::coords(rocobj, "best"))$threshold
  trt_assign2$dopt <- with(trt_assign2, ifelse(prob_dopt <= yindex2, "dopt=0",
                                               "dopt=1"))
  
  # Run SUBMOD #
  mod1 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="mob_weib", pool = "otr:logistic", 
                       hyper = list(maxdepth=3), delta = delta,
                       otr_prob_thres = "youden")
  
  eql_otrlog <- all.equal(trt_assign1[,c("Subgrps", "dopt")], 
                          mod1$trt_assign)
  
  mod2 <- submod_train(Y=Y, A=A, X=X, family = family,
                       mu_train = mu_train,
                       submod="mob_weib", pool = "otr:rf", 
                       hyper = list(maxdepth=3), delta = delta,
                       otr_prob_thres = "youden")
  
  eql_otrrf <- all.equal(trt_assign2[,c("Subgrps", "dopt")], 
                         mod2$trt_assign)
  
  # Check if the youden index's are similar #
  diff_yindex1 <- abs(yindex1 - unique(mod1$trt_eff_pool$thres_youden))
  diff_yindex2 <- abs(yindex2 - unique(mod2$trt_eff_pool$thres_youden))
  
  bound <- 0.05
  diff_yindex1 <- ifelse(diff_yindex1 < bound, TRUE, FALSE)
  diff_yindex2 <- ifelse(diff_yindex2 < bound, TRUE, FALSE)
  
  # Check #
  expect_equal(eql_trteff, TRUE)
  expect_equal(eql_otrlog, TRUE)
  expect_equal(eql_otrrf, TRUE)
  expect_equal(diff_yindex1, TRUE)
  expect_equal(diff_yindex2, TRUE)
  
})
