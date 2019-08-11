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
#'   \item mods - trained model(s)
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
#' summary(mod1$mu_train$PLE)
#'

##### Elastic net (glmnet): Y~(A,X,A*X) ==> PLEs ######
ple_glmnet = function(Y, A, X, Xtest, lambda="lambda.min", family, ...){

  ## Model matrix the covariate space ###
  X = model.matrix(~. -1, data = X )
  ## Generate the interaction covariates ###
  X_inter = X*A
  colnames(X_inter) = paste(colnames(X), "_A", sep="")
  W = cbind(X, A, X_inter)
  ##### Elastic Net #####
  if (family=="survival") { family = "cox"  }
  mod <- cv.glmnet(x = W, y = Y, alpha=0.5, family=family)
  res = list(mods = mod, lambda=lambda, X=X)
  class(res) = "ple_glmnet"
  ## Return Results ##
  return( res )
}

#' Predict Patient-level Estimates: glmnet
#'
#' For continuous/binary (family="gaussian" or "binomial"), output estimates of
#' (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0)). For survival, output estimates of
#' (HR(X,A=1), HR(X,A=0), HR(X, A=1)-HR(X, A=0)).
#'
#' @param object Trained glmnet model(s).
#' @param newdata Data-set to make predictions at (Default=NULL, predictions correspond
#' to training data).
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import glmnet
#'
#' @return Data-frame with predictions of (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0))
#' for continuous/binary outcomes. For survival, returns (HR(X,A=1), HR(X,A=0),
#' HR(X, A=1)-HR(X, A=0)).
#'
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
#' summary(mod1$mu_train)
#' summary(predict(mod1)) # Training set predictions, newdata=NULL
#' summary(predict(mod1, newdata=X)) # Test set predictions
#'
#' @method predict ple_glmnet
#' @export
#'
#### Predict: ple_glmnet ####
predict.ple_glmnet = function(object, newdata=NULL, ...){

  mod = object$mods
  lambda = object$lambda
  X = object$X
  if (!is.null(newdata)){
    X = newdata
    X = model.matrix(~. -1, data = X )
  }
  ### Predictions (Counterfactuals) ###
  mu_hat = data.frame( mu1 =  as.numeric(predict(mod,
                                                   newx = cbind(X, 1, X*1), s=lambda)),
                       mu0 =  as.numeric(predict(mod,
                                                   newx = cbind(X, 0, X*0), s=lambda)) )
  mu_hat$PLE = with(mu_hat, mu1 - mu0 )
  return( mu_hat  )
}
