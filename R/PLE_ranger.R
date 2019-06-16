#' Patient-level Estimates: Ranger
#'
#' Estimate patient-level estimates (PLEs) through treatment-specific
#' (or with explicit X*A interactions) random forest models (ranger).
#' Used for continuous or binary outcomes, with output estimates of E(Y|X,A=a) and
#' E(Y|X,A=1)-E(Y|X,A=0) (PLE).
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param byTrt If 1, fit treatment-specific ranger models. If 0, fit a single ranger model
#' with covariate space (X, A, X*A)
#' @param min.node.pct Minimum sample size in forest nodes (n*min.node.pct)
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import ranger
#'
#' @return Patient-level estimates (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0)) for train/test sets.
#'  \itemize{
#'   \item mods - trained model(s)
#'   \item mu_train - Patient-level estimates (training set)
#'   \item mu_test - Patient-level estimates (test set)
#' }
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' \donttest{
#' # Default (treatment-specific ranger models) #
#' mod1 = PLE_ranger(Y, A, X, Xtest=X)
#' summary(mod1$mu_train$PLE)
#'
#' # Generate A*X covariates (single ranger model) #
#' mod2 = PLE_ranger(Y, A, X, Xtest=X, byTrt=0)
#' summary(mod2$mu_train$PLE)
#' }
#'
#' ## Survival (TBD) ##
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{ranger}}
#'
#### Counterfactual Forest: Ranger ####
PLE_ranger = function(Y, A, X, Xtest, byTrt=1, min.node.pct=0.10, family="gaussian", ...){

  set.seed(668877)
  ## Random Forest models for each Treatment ##
  if (byTrt==1){
    ## Split data by treatment ###
    train0 =  data.frame(Y=Y[A==0], X[A==0,])
    train1 =  data.frame(Y=Y[A==1], X[A==1,])
    # Trt 0 #
    mod0 <- ranger(Y ~ ., data = train0, seed=1, min.node.size = min.node.pct*dim(train0)[1])
    # Trt 1 #
    mod1 <- ranger(Y ~ ., data = train1, seed=2, min.node.size = min.node.pct*dim(train1)[1])

    mods = list(mod0=mod0, mod1=mod1)
    if (family %in% c("gaussian", "binomial")){
      ## Predictions: Train/Test ##
      mu_train = data.frame( mu1 =  predict(mod1, data = X)$predictions,
                             mu0 = predict(mod0, data = X)$predictions)
      mu_train$PLE = with(mu_train, mu1 - mu0 )

      mu_test = data.frame( mu1 =  predict(mod1, data = Xtest)$predictions,
                            mu0 = predict(mod0, data = Xtest)$predictions)
      mu_test$PLE = with(mu_test, mu1 - mu0 )
    }
    if (family=="survival"){
      ## Predict Difference in RMST (can restrict up to time-point) ##
      ## TBD ##
    }
  }
  ## Single Random Forest Model: Generate A*X interactions manually ##
  if (byTrt==0){

    ## Set up A*X interactions in training set ##
    X_inter = X*A
    colnames(X_inter) = paste(colnames(X), "_A", sep="")
    train.inter = cbind(Y, A, X, X_inter)
    ## Set up Counter-factual data-sets for train/test ##
    Xtrain0 = data.frame(0, X, X*0)
    colnames(Xtrain0) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
    Xtrain1 = data.frame(1, X, X*1)
    colnames(Xtrain1) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
    ## A*X interactions in test set; for A=1 and A=0 separately ##
    Xtest0 = data.frame(0, Xtest, Xtest*0)
    colnames(Xtest0) = c( "A", colnames(Xtest), paste(colnames(Xtest), "_A", sep="") )
    Xtest1 = data.frame(1, Xtest, Xtest*1)
    colnames(Xtest1) = c( "A", colnames(Xtest), paste(colnames(Xtest), "_A", sep="") )

    ## Fit RF ##
    mod.inter <- ranger(Y ~ ., data = train.inter, seed=5,
                 min.node.size = min.node.pct*dim(train.inter)[1])

    ## Predictions: Train/Test ##
    mu_train = data.frame( mu1 =  predict(mod.inter, data = Xtrain1)$predictions,
                           mu0 = predict(mod.inter, data = Xtrain0)$predictions)
    mu_train$PLE = with(mu_train, mu1 - mu0 )

    mu_test = data.frame( mu1 =  predict(mod.inter, data = Xtest1)$predictions,
                          mu0 = predict(mod.inter, data = Xtest0)$predictions)
    mu_test$PLE = with(mu_test, mu1 - mu0 )
    mods = list(mod.inter=mod.inter)
  }

  ## Return Results ##
  return( list(mods=mods, mu_train = mu_train, mu_test = mu_test) )
}
