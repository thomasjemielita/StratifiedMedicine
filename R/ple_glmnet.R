#' Patient-level Estimates: Elastic Net (glmnet)
#'
#' Uses the elastic net (glmnet R package) to obtain patient-level estimates. Usable for
#' continuous, binary, or survival outcomes.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate space.
#' @param Xtest Test set
#' @param lambda Lambda for elastic-net (default = "lambda.min"). Other options include
#' "lambda.1se" or  fixed values
#' @param family Outcome type ("gaussian", "binomial", "survival"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import glmnet
#'
#' @return Trained glmnet model(s).
#'  \itemize{
#'   \item mod - trained model(s)
#'   \item lambda - Lambda used for elastic-net (passes to prediction function)
#'   \item X - Covariate Space (in model matrix form)
#' }
#' @export
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#'
#' mod1 = ple_glmnet(Y, A, X, Xtest=X, family="gaussian")
#'

##### Elastic net (glmnet): Y~(A,X,A*X) ==> PLEs ######
ple_glmnet = function(Y, A, X, Xtest, lambda="lambda.min", family, ...){

  ## Model matrix the covariate space ###
  X = model.matrix(~. -1, data = X )
  if (is.null(A)){
    W = X
  }
  if (!is.null(A)){
    ## Generate the interaction covariates ###
    X_inter = X*A
    colnames(X_inter) = paste(colnames(X), "_A", sep="")
    W = cbind(X, A, X_inter)
  }
  ##### Elastic Net #####
  if (family=="survival") { family = "cox"  }
  mod <- cv.glmnet(x = W, y = Y, alpha=0.5, family=family)
  mod <- list(mod=mod, lambda=lambda, A = length(A) )
  pred.fun = function(mod, X=NULL){
    lambda = mod$lambda
    A = mod$A
    mod = mod$mod
    if (!is.null(X)){
      X = model.matrix(~. -1, data = X )
    }
    ### Predictions (Counterfactuals) ###
    if (A==0){
      mu_hat = data.frame( mu=as.numeric(predict(mod, newx = X, s=lambda)) ) 
    }
    if (A>0){
      mu_hat = data.frame( 
        mu1 = as.numeric(predict(mod,newx = cbind(X, 1, X*1), s=lambda)),
        mu0 = as.numeric(predict(mod,newx = cbind(X, 0, X*0), s=lambda)) )
      mu_hat$PLE = with(mu_hat, mu1 - mu0 )
    }
    return(mu_hat)
  }
  res = list(mod = mod, pred.fun=pred.fun)
  class(res) = "ple_glmnet"
  ## Return Results ##
  return( res )
}