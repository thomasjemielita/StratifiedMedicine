test_that("Test whether submod_train works (ctns)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  
  dat_ctns = generate_subgrp_data(family="gaussian")
  Y = dat_ctns$Y
  X = dat_ctns$X
  A = dat_ctns$A
  family <- "gaussian"
  
  # lmtree #
  mod1 <- submod_train(Y=Y, A=A, X=X, family = family, submod="lmtree")
  
  sub1 <- lmtree(Y~A | . , data = X, alpha = 0.10, 
                 maxdepth = 4, minsize = floor(dim(X)[1]*0.10))
  
  pred_s1 <- as.character(predict(sub1, type="node", newdata = X))
  
  eql_loglik <- all.equal(logLik(mod1$fit$mod), logLik(sub1))
  eql_mob <- all.equal(mod1$Subgrps.train, pred_s1)
  
  # Counterfactual based #
  ple_mod <- ple_train(Y=Y, A=A, X=X, family = family)
  mu_train <- ple_mod$mu_train
  
  # CTREE (counterfactuals) #
  mod2 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="ctree_cate")
  
  sub2 <- ctree(mu_train$diff_1_0 ~ ., data = X,
                control = partykit::ctree_control(alpha=0.10, 
                                                  minbucket=floor(dim(X)[1]*0.10), 
                                                  maxdepth=4))
  pred_s2 <- as.character(predict(sub2, type="node", newdata = X))
  
  eql_ctree <- all.equal(mod2$Subgrps.train, pred_s2)
  
  # RPART (counterfactuals) #
  mod3 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="rpart_cate")
  
  sub3 <- rpart(mu_train$diff_1_0 ~ ., data = X,
                control = rpart::rpart.control(minbucket=floor(dim(X)[1]*0.10), 
                                                  maxdepth=4))
  sub3 <- as.party(sub3)
  pred_s3 <- as.character(predict(sub3, type="node", newdata = X))
  
  eql_rpart <- all.equal(mod3$Subgrps.train, pred_s3)
  
  
  # OTR (weighted ctree) #
  mod4 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="otr", delta = "> 0.10")
  
  ind_PLE <- ifelse(mu_train$diff_1_0 > 0.10, 1, 0)
  w_PLE <- abs(mu_train$diff_1_0 - 0.10)
  hold <- data.frame(ind_PLE, X)
  sub4 <- ctree(ind_PLE ~ ., data = hold, weights = w_PLE, 
               control = partykit::ctree_control(alpha=0.10, 
                                                 minbucket=floor(dim(X)[1]*0.10),
                                                 maxdepth=4))
  pred_s4 <-as.character(predict(sub4, type="node", newdata = X))

  eql_otr <- all.equal(mod4$Subgrps.train, pred_s4)
  
  
  # Check #
  expect_equal(eql_loglik, TRUE)
  expect_equal(eql_mob, TRUE)
  expect_equal(eql_ctree, TRUE)
  expect_equal(eql_rpart, TRUE)
  expect_equal(eql_otr, TRUE)

})

test_that("Test whether submod_train works (binomial)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  
  dat_bin = generate_subgrp_data(family="binomial")
  Y = dat_bin$Y
  X = dat_bin$X
  A = dat_bin$A
  family <- "binomial"
  
  # glmtree #
  mod1 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       submod="glmtree", hyper = list(maxdepth=3))
  
  sub1 <- glmtree(Y~A | . , data = X, alpha = 0.10, 
                  family= binomial(link="identity"), 
                 maxdepth = 3, minsize = floor(dim(X)[1]*0.10))
  
  pred_s1 <- as.character(predict(sub1, type="node", newdata = X))
  
  eql_loglik <- all.equal(logLik(mod1$fit$mod), logLik(sub1))
  eql_mob <- all.equal(mod1$Subgrps.train, pred_s1)
  
  # Counterfactual based #
  ple_mod <- ple_train(Y=Y, A=A, X=X, family = family)
  mu_train <- ple_mod$mu_train
  
  # CTREE (counterfactuals) #
  mod2 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="ctree_cate", 
                       hyper = list(maxdepth=3))
  
  sub2 <- ctree(mu_train$diff_1_0 ~ ., data = X,
                control = partykit::ctree_control(alpha=0.10, 
                                                  minbucket=floor(dim(X)[1]*0.10), 
                                                  maxdepth=3))
  pred_s2 <- as.character(predict(sub2, type="node", newdata = X))
  
  eql_ctree <- all.equal(mod2$Subgrps.train, pred_s2)
  
  # RPART (counterfactuals) #
  mod3 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="rpart_cate", 
                       hyper = list(maxdepth=3))
  
  sub3 <- rpart(mu_train$diff_1_0 ~ ., data = X,
                control = rpart::rpart.control(minbucket=floor(dim(X)[1]*0.10), 
                                               maxdepth=3))
  sub3 <- as.party(sub3)
  pred_s3 <- as.character(predict(sub3, type="node", newdata = X))
  
  eql_rpart <- all.equal(mod3$Subgrps.train, pred_s3)
  
  
  # OTR (weighted ctree) #
  mod4 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="otr", delta = "> 0.02",
                       hyper = list(maxdepth=3))
  
  ind_PLE <- ifelse(mu_train$diff_1_0 > 0.02, 1, 0)
  w_PLE <- abs(mu_train$diff_1_0 - 0.02)
  hold <- data.frame(ind_PLE, X)
  sub4 <- ctree(ind_PLE ~ ., data = hold, weights = w_PLE, 
                control = partykit::ctree_control(alpha=0.10, 
                                                  minbucket=floor(dim(X)[1]*0.10),
                                                  maxdepth=3))
  pred_s4 <- as.character(predict(sub4, type="node", newdata = X))
  
  eql_otr <- all.equal(mod4$Subgrps.train, pred_s4)
  
  
  # Check #
  expect_equal(eql_mob, TRUE)
  expect_equal(eql_loglik, TRUE)
  expect_equal(eql_ctree, TRUE)
  expect_equal(eql_rpart, TRUE)
  expect_equal(eql_otr, TRUE)
  
})

test_that("Test whether submod_train works (survival)", {
  
  skip_on_cran()
  
  library(partykit)
  library(rpart)
  library(survival)
  require(TH.data); require(coin)
  data("GBSG2", package = "TH.data")
  surv.dat = GBSG2
  # Design Matrices ###
  Y = with(surv.dat, Surv(time, cens))
  X = surv.dat[,!(colnames(surv.dat) %in% c("time", "cens")) ]
  set.seed(513)
  A = rbinom(n = dim(X)[1], size=1, prob=0.5)
  
  family <- "survival"
  
  # LMTREE #
  mod0 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       submod="lmtree", hyper = list(maxdepth=3))
  
  Y_use <- coin::logrank_trafo(Y)
  sub0 <- partykit::lmtree(Y_use ~ A | ., data = X, 
                           alpha=0.10, minsize=floor(dim(X)[1]*0.10), 
                           maxdepth=3)
  
  pred_s0 <- as.character(predict(sub0, type="node", newdata = X))
  
  eql_loglik_lmtree <- all.equal(logLik(mod0$fit$mod), logLik(sub0))
  eql_lmtree <- all.equal(mod0$Subgrps.train, pred_s0)
  
  # MOB-WEIB #
  mod1 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       submod="mob_weib", hyper = list(maxdepth=3))
  
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
  
  pred_s1 <- as.character(predict(sub1, type="node", newdata = X))
  
  eql_loglik <- all.equal(logLik(mod1$fit$mod), logLik(sub1))
  eql_mob <- all.equal(mod1$Subgrps.train, pred_s1)
  
  # Counterfactual based #
  ple_mod <- ple_train(Y=Y, A=A, X=X, family = family)
  mu_train <- ple_mod$mu_train
  
  # CTREE (counterfactuals) #
  mod2 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="ctree_cate", 
                       hyper = list(maxdepth=3))
  
  sub2 <- ctree(mu_train$diff_1_0 ~ ., data = X,
                control = partykit::ctree_control(alpha=0.10, 
                                                  minbucket=floor(dim(X)[1]*0.10), 
                                                  maxdepth=3))
  pred_s2 <-as.character(predict(sub2, type="node", newdata = X))
  
  eql_ctree <- all.equal(mod2$Subgrps.train, pred_s2)
  
  # RPART (counterfactuals) #
  mod3 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="rpart_cate", 
                       hyper = list(maxdepth=3))
  
  sub3 <- rpart(mu_train$diff_1_0 ~ ., data = X,
                control = rpart::rpart.control(minbucket=floor(dim(X)[1]*0.10), 
                                               maxdepth=3))
  sub3 <- as.party(sub3)
  pred_s3 <- as.character(predict(sub3, type="node", newdata = X))
  
  eql_rpart <- all.equal(mod3$Subgrps.train, pred_s3)
  
  
  # OTR (weighted ctree) #
  mod4 <- submod_train(Y=Y, A=A, X=X, family = family, 
                       mu_train = mu_train,
                       submod="otr", delta = "> 50",
                       hyper = list(maxdepth=3))
  
  ind_PLE <- ifelse(mu_train$diff_1_0 > 50, 1, 0)
  w_PLE <- abs(mu_train$diff_1_0 - 50)
  hold <- data.frame(ind_PLE, X)
  sub4 <- ctree(ind_PLE ~ ., data = hold, weights = w_PLE, 
                control = partykit::ctree_control(alpha=0.10, 
                                                  minbucket=floor(dim(X)[1]*0.10),
                                                  maxdepth=3))
  pred_s4 <- as.character(predict(sub4, type="node", newdata = X))
  
  eql_otr <- all.equal(mod4$Subgrps.train, pred_s4)
  
  
  # Check #
  expect_equal(eql_lmtree, TRUE)
  expect_equal(eql_loglik_lmtree, TRUE)
  expect_equal(eql_mob, TRUE)
  expect_equal(eql_loglik, TRUE)
  expect_equal(eql_ctree, TRUE)
  expect_equal(eql_rpart, TRUE)
  expect_equal(eql_otr, TRUE)
  
})
