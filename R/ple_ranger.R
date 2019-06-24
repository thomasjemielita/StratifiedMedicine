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
#' @return Trained random forest (ranger) model(s).
#'  \itemize{
#'   \item mods - trained model(s)
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
#' mod1 = ple_ranger(Y, A, X, Xtest=X)
#' summary( predict(mod1, newdata=data.frame(A,X) ) ) # oob predictions for training
#' summary( predict(mod1, newdata=data.frame(X) ) ) # new-predictions, no oob here
#'
#' # Generate A*X covariates (single ranger model) #
#' mod2 = ple_ranger(Y, A, X, Xtest=X, byTrt=0)
#' summary( predict(mod2, newdata=data.frame(A,X) ) ) # oob predictions for training
#' summary( predict(mod2, newdata=data.frame(X) ) ) # new-predictions
#' }
#'
#' ## Survival (TBD) ##
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{ranger}}
#'
#### Counterfactual Forest: Ranger ####
ple_ranger = function(Y, A, X, Xtest, byTrt=1, min.node.pct=0.10, family="gaussian", ...){

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
    mods = list(mod.inter=mod.inter)
  }
  res = list(mods=mods)
  class(res) = "ple_ranger"
  ## Return Results ##
  return( res )
}

#' Predict Patient-level Estimates: Ranger
#'
#' Get estimates of (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0)) using trained random
#' forest (ranger) model(s).
#'
#' @param object Trained random forest (ranger) model(s).
#' @param newdata Data-set to make predictions at. For training data, should
#' include covariates (X) and treatment (A) for oob predictions.
#' @param oob Use out-of-bag predictions (default=TRUE unless newdata does not contain
#' treatment variable A).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import ranger
#'
#' @return Data-frame with predictions of (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0))
#'
#' @examples
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#'
#' # Default (treatment-specific ranger models) #
#' mod1 = ple_ranger(Y, A, X, Xtest=X)
#' summary( predict(mod1, newdata=data.frame(A,X) ) ) # oob predictions for training
#' summary( predict(mod1, newdata=data.frame(X) ) ) # new-predictions, no oob here
#'
#'
#' @method predict ple_ranger
#' @export
#'
#### Counterfactual Forest: Ranger ####
predict.ple_ranger = function(object, newdata, oob=TRUE, ...){

  if ( !("A" %in% names(newdata)) ){
    oob = FALSE
  }
  A = newdata$A
  X = newdata[,!(colnames(newdata) %in% "A") ]
  ## Treatment-specific ranger models ##
  if (length(object$mods)>1){
    mu1_hat = predict( object$mods[["mod1"]], X )$predictions
    mu0_hat = predict( object$mods[["mod0"]], X )$predictions
    if (oob){
      mu1_hat = ifelse(A==1, object$mods[["mod1"]]$predictions, mu1_hat)
      mu0_hat = ifelse(A==0, object$mods[["mod0"]]$predictions, mu0_hat)
    }
  }
  ## X*A interactions included (one ranger model) ##
  if (length(object$mods)==1){
    X0 = data.frame(0, X, X*0)
    colnames(X0) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
    X1 = data.frame(1, X, X*1)
    colnames(X1) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
    mu1_hat = predict( object$mods[["mod.inter"]], X1 )$predictions
    mu0_hat = predict( object$mods[["mod.inter"]], X0 )$predictions
    if (oob){
      mu1_hat = ifelse(A==1, object$mods[["mod.inter"]]$predictions, mu1_hat)
      mu0_hat = ifelse(A==0, object$mods[["mod.inter"]]$predictions, mu0_hat)
    }
  }
  mu_hat = data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = mu1_hat-mu0_hat)
  return( mu_hat  )
}

