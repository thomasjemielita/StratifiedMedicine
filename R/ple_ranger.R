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
#' model with covariate space (X, A, X*A). For "gaussian" or "binomial", default is TRUE.
#' For "survival", default is FALSE. 
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
#'
#' }
#'
#' @export
#' @seealso \code{\link{PRISM}}, \code{\link{ranger}}
#'
#### Counterfactual Forest: Ranger ####
ple_ranger = function(Y, A, X, Xtest, family="gaussian", 
                      byTrt=FALSE,
                      min.node.pct=0.10,...){

  if (is.Surv(Y) & family!="survival"){
    byTrt = FALSE
  }
  if (is.null(A)){
    mod <- ranger(Y ~ ., data = data.frame(Y, X), 
                   min.node.size = min.node.pct*dim(X)[2] )
    mod = list(mod=mod)
    # Prediction Function #
    pred.fun = function(mod, X, tau=NULL){
      treetype = mod[[1]]$treetype
      if (treetype!="Survival"){
          mu_hat = data.frame(PLE = predict( mod$mod, X )$predictions )
      }
      if (treetype=="Survival"){
        preds = predict( mod$mod, X )
        if (is.null(tau)){
          tau.t <- max(preds$unique.death.times)
        }
        looper_rmst <- function(i, surv, time){
          est.rmst <- rmst_calc(surv = surv[i,],
                                time = time,
                                tau=tau.t)
          return(est.rmst)
        }
        rmst <- lapply(1:dim(X)[1], looper_rmst, surv=preds$survival,
                       time=preds$unique.death.times)
        mu_hat <- do.call(rbind, rmst)
        mu_hat <- data.frame(PLE = mu_hat)
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
      mod0 <- ranger(Y ~ ., data = train0, min.node.size = min.node.pct*dim(train0)[1])
      # Trt 1 #
      mod1 <- ranger(Y ~ ., data = train1, min.node.size = min.node.pct*dim(train1)[1])
      mod = list(mod0=mod0, mod1=mod1)
      # Prediction Function #
      pred.fun = function(mod, X, tau=NULL){
        treetype = mod[[1]]$treetype
        if (treetype!="Survival"){
          mu1_hat = predict( mod$mod1, X )$predictions
          mu0_hat = predict( mod$mod0, X )$predictions
          mu_hat = data.frame(mu1 = mu1_hat, mu0 = mu0_hat, PLE = mu1_hat-mu0_hat) 
        }
        if (treetype=="Survival"){
          pred1 = predict( mod$mod1, X ) 
          pred0 = predict( mod$mod0, X )
          if (is.null(tau)){
            tau.t <- min( max(pred0$unique.death.times), max(pred1$unique.death.times))
          }
          looper_rmst <- function(i, surv, time){
            est.rmst <- rmst_calc(surv = surv[i,],
                                  time = time,
                                  tau=tau.t)
            return(est.rmst)
          }
          rmst1 <- lapply(1:dim(X)[1], looper_rmst, surv=pred1$survival,
                         time=pred1$unique.death.times)
          rmst1 <- do.call(rbind, rmst1)
          rmst0 <- lapply(1:dim(X)[1], looper_rmst, surv=pred0$survival,
                          time=pred0$unique.death.times)
          rmst0 <- do.call(rbind, rmst0)
          mu_hat <- data.frame(mu1 = rmst1, mu0 = rmst0, PLE=rmst1-rmst0)
        }
        return(mu_hat)
      }
    }
    ## Single Random Forest Model: Generate A*X interactions manually ##
    if (!byTrt){
      X = model.matrix(~., data = X )
      X = X[,-1]
      ## Set up A*X interactions in training set ##
      X_inter = X*A
      colnames(X_inter) = paste(colnames(X), "_A", sep="")
      train.inter = data.frame(Y, A, X, X_inter)
      ## Fit RF ##
      mod.inter <- ranger(Y ~ ., data = train.inter,
                          min.node.size = min.node.pct*dim(train.inter)[1])
      mod = list(mod.inter=mod.inter)
      # Prediction Function #
      pred.fun = function(mod, X, tau=NULL){
        X = model.matrix(~., data = X )
        X = X[,-1]
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
          pred1 = predict( mod$mod.inter, X1 ) 
          pred0 = predict( mod$mod.inter, X0 )
          if (is.null(tau)){
            tau.t <- min( max(pred0$unique.death.times), 
                          max(pred1$unique.death.times)  )
          }
          looper_rmst <- function(i, surv, time){
            est.rmst <- rmst_calc(surv = surv[i,],
                                  time = time,
                                  tau=tau.t)
            return(est.rmst)
          }
          # A = 0 #
          rmst0 <- lapply(1:dim(X)[1], looper_rmst, surv=pred0$survival,
                          time=pred0$unique.death.times)
          rmst0 <- do.call(rbind, rmst0)
          # A = 1 #
          rmst1 <- lapply(1:dim(X)[1], looper_rmst, surv=pred1$survival,
                          time=pred1$unique.death.times)
          rmst1 <- do.call(rbind, rmst1)
          mu_hat <- data.frame(mu1 = rmst1, mu0 = rmst0)
          mu_hat$PLE <- with(mu_hat, mu1 - mu0)
        }
        return(mu_hat)
      }
    }
  }
  res = list(mod=mod, pred.fun=pred.fun, A=A, X=X)
  class(res) = "ple_ranger"
  ## Return Results ##
  return( res )
}