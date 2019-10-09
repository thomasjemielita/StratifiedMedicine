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
#'   \item mod - trained model(s)
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
    mod = list(mod=mod)
    # Prediction Function #
    pred.fun = function(mod, X){
      treetype = mod[[1]]$treetype
      if (treetype!="Survival"){
          mu_hat = data.frame(PLE = predict( mod$mod, X )$predictions )
      }
      if (treetype=="Survival"){
        preds = predict( mod$mod, X )
        mu_hat = ranger_rmst(preds, X, trt=FALSE)
      }
      return(mu_hat)
    }
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
      mod = list(mod0=mod0, mod1=mod1)
      # Prediction Function #
      pred.fun = function(mod, X){
        treetype = mod[[1]]$treetype
        if (treetype!="Survival"){
          mu1_hat = predict( mod$mod1, X )$predictions
          mu0_hat = predict( mod$mod0, X )$predictions
          mu_hat = data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = mu1_hat-mu0_hat) 
        }
        if (treetype=="Survival"){
          pred1 = predict( mod$mod1, X ) 
          pred0 = predict( mod$mod0, X )
          mu_hat = ranger_rmst(preds=list(pred1=pred1,pred0=pred0), X, trt=TRUE)
        }
        return(mu_hat)
      }
    }
    ## Single Random Forest Model: Generate A*X interactions manually ##
    if (!byTrt){
      ## Set up A*X interactions in training set ##
      X_inter = X*A
      colnames(X_inter) = paste(colnames(X), "_A", sep="")
      train.inter = data.frame(Y, A, X, X_inter)
      ## Fit RF ##
      mod.inter <- ranger(Y ~ ., data = train.inter, seed=5,
                          min.node.size = min.node.pct*dim(train.inter)[1])
      mod = list(mod.inter=mod.inter)
      # Prediction Function #
      pred.fun = function(mod, X){
        mod.inter = mod$mod.inter
        treetype = mod$mod.inter$treetype
        X0 = data.frame(0, X, X*0)
        colnames(X0) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
        X1 = data.frame(1, X, X*1)
        colnames(X1) = c( "A", colnames(X), paste(colnames(X), "_A", sep="") )
        if (treetype!="Survival"){
          mu_hat = data.frame(mu1 = predict(mod$mod.inter, X1)$predictions,
                              mu0 = predict(mod$mod.inter, X0)$predictions )
          mu_hat$PLE = with(mu_hat, mu1-mu0)
        }
        if (treetype=="Survival"){
          pred1 = predict( mod$mod.inter, data=X1) 
          pred0 = predict( mod$mod.inter, data=X0)
          mu_hat = ranger_rmst(preds=list(pred1=pred1, pred0=pred0), X, trt=TRUE)
        }
        return(mu_hat)
      }
      summary(pred.fun(mod, X))
    }
  }
  res = list(mod=mod, pred.fun=pred.fun, A=A, X=X)
  class(res) = "ple_ranger"
  ## Return Results ##
  return( res )
}