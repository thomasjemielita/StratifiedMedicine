#' Patient-level Estimates: Ranger
#'
#' Uses treatment-specific (or with explicit X*A interactions) random forest models (ranger)
#' to obtain patient-level estimates. Used for continuous, binary, or survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param byTrt If TRUE, fit treatment-specific ranger models. If FALSE, fit a single ranger
#' model with covariate space (X, A, X*A).
#' @param min.node.pct Minimum sample size in forest nodes (n*min.node.pct)
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import ranger
#'
#' @return Trained random forest (ranger) model(s).
#'  \itemize{
#'   \item mods - trained model(s)
#'   \item A - treatment variable (training set)
#'   \item X - covariate space (training set)
#' }
#' @examples
#' \donttest{
#' library(StratifiedMedicine)
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' # Counter-factual Random Forest (treatment-specific ranger models) #
#' mod1 = ple_ranger(Y, A, X, Xtest=X)
#' summary( predict(mod1, newdata=data.frame(A,X) ) ) # oob predictions for training
#' summary( predict(mod1, newdata=data.frame(X) ) ) # new-predictions, no oob here
#'
#' }
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{ranger}}
#'
#### Counterfactual Forest: Ranger ####
ple_ranger = function(Y, A, X, Xtest, byTrt=TRUE, min.node.pct=0.10, family="gaussian",
                      ...){

  if (is.null(A)){
    mod <- ranger(Y ~ ., data = data.frame(Y, X), seed=1, 
                   min.node.size = min.node.pct*dim(X)[2] )
    mods = list(mod=mod)
  }
  if (!is.null(A)){
    ## Random Forest models for each Treatment ##
    if (byTrt){
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
    if (!byTrt){
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
  }
  res = list(mods=mods, A=A, X=X)
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
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param oob Use out-of-bag predictions (default=TRUE). Only applicable for training data
#' (newdata=NULL).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import ranger
#'
#' @return Data-frame with predictions of (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0)) or
#' survival probabilities and difference in restricted mean survival time (RMST),
#' (S(T|X,A=1), S(T|X,A=0), RMST(A=1,X)-RMST(A=0,X) )
#'
#' @examples
#' \donttest{
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
#' summary( predict(mod1 ) ) # oob predictions for training
#' summary( predict(mod1, newdata=X ) ) # new-predictions, no oob here
#' }
#'
#' @method predict ple_ranger
#' @export
#'
#### Counterfactual Forest: Ranger ####
predict.ple_ranger = function(object, newdata=NULL, oob=FALSE, ...){

  mods = object$mods
  A = object$A
  X = object$X
  if (!is.null(newdata)){
    oob=FALSE
    X = newdata
  }
  treetype = mods[[1]]$treetype
  if (is.null(A)){
    if (treetype!="Survival"){
      mu_hat = data.frame(PLE = predict( object$mods$mod, X )$predictions )
    } 
    if (treetype=="Survival"){
      preds = predict( mods$mod, X )
      mu_hat = data.frame( preds$survival )
      colnames(mu_hat) = paste("time", mods$mod$unique.death.times)
    }
  }
  if (!is.null(A)){
    ### E(Y|X,A=1) - E(Y|X,A=0) ###
    if (treetype!="Survival"){
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
      ple_hat = mu1_hat-mu0_hat
    }
    ## Predict Difference in RMST (can restrict up to time-point) ##
    if (treetype=="Survival"){
      if (!requireNamespace("zoo", quietly = TRUE)) {
        stop("Package zoo needed for ranger RMST predictions. Please install.")
      }
      if (length(object$mods)>1){
        out0 = predict(mods[["mod0"]], data=X)
        out1 = predict(mods[["mod1"]], data=X)
      }
      if (length(object$mods)==1){
        X0 = data.frame(0, X, X*0)
        colnames(X0) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
        X1 = data.frame(1, X, X*1)
        colnames(X1) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
        out0 = predict(mods[["mod.inter"]], data=X)
        out1 = predict(mods[["mod.inter"]], data=X)
      }
      times0 = out0$unique.death.times
      mu0_hat = out0$survival
      times1 = out1$unique.death.times
      mu1_hat = out1$survival
      # Order times and calculate RMST (A=1 vs A=0) #
      id1 <- order(times1)
      id0 <- order(times0)
      ple_hat = NULL
      dim.length = ifelse(is.null(newdata), dim(X)[1], dim(newdata)[1])
      for (ii in 1:dim.length){
        ## Trt 1 ##
        y <- mu1_hat[ii,]
        rmst1 <- sum(diff(times1[id1])*zoo::rollmean(y[id1],2))
        ## Trt 0 ##
        y <- mu0_hat[ii,]
        rmst0 <- sum(diff(times0[id0])*zoo::rollmean(y[id0],2))
        # Store #
        ple_hat = c(ple_hat, (rmst1-rmst0) )
      }
      ## take average of individual survival probabilities for mu1_hat and mu0_hat ##
      mu1_hat = apply(mu1_hat, 1, mean)
      mu0_hat = apply(mu0_hat, 1, mean)
    }
    mu_hat = data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = ple_hat)
  }

  return( mu_hat  )
}

