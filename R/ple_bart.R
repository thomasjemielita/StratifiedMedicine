#' Patient-level Estimates: BART
#'
#' Uses the BART algorithm (from BART R package) to obtain patient-level estimates.
#' Used for continuous or binary outcomes, with output estimates of
#' E(Y|X,A=a) and E(Y|X,A=1)-E(Y|X,A=0) (PLE). In progress.
#'
#' @param Y The outcome variable. Must be numeric or survival (ex; Surv(time,cens) )
#' @param A Treatment variable. (a=1,...A)
#' @param X Covariate matrix. Must be numeric.
#' @param Xtest Test set
#' @param family Outcome type ("gaussian", "binomial"), default is "gaussian"
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import BART
#'
#' @return Trained BART model(s) and patient-level estimates
#' (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0)) for train/test sets.
#'  \itemize{
#'   \item mods - trained model(s)
#'   \item mu_train - Patient-level estimates (training set)
#'   \item mu_test - Patient-level estimates (test set)
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
#' train = data.frame(Y, A, X)
#'
#' # BART #
#' \donttest{
#' mod1 = ple_bart(train=train, Xtest=X)
#' summary(mod1$mu_train)
#' summary(predict(mod1, newdata=X))
#'
#' summary(mod1$mu_train$PLE)
#' }
#'
#' @seealso \code{\link{PRISM}}, \code{\link{BART}}
#'
#### BART ###
ple_bart = function(Y, A, X, Xtest, family="gaussian", ...){

  ## Generate counterfactual design matrices  ##
  Xtrain_0 = data.frame(A=0, X)
  Xtrain_1 = data.frame(A=1, X)
  Xtest_0 = data.frame(A=0, Xtest)
  Xtest_1 = data.frame(A=1, Xtest)

  x.train = data.frame(A,X)
  # Test: Include both train/test for faster predictions #
  X.FULL = rbind(Xtrain_0, Xtrain_1, Xtest_0, Xtest_1)
  ## BART ##
  set.seed(51351)
  bartFit = wbart(x.train = x.train, y.train = Y, x.test = X.FULL)
  n.tr = dim(Xtrain_0)[1]
  n.ts = dim(Xtest_0)[1]
  ### PLE Predictions: Train/Test ###
  mu_train = data.frame(mu1 = bartFit$yhat.test.mean[(n.tr+1):(2*n.tr)],
                        mu0 = bartFit$yhat.test.mean[1:n.tr])
  mu_train$PLE = mu_train$mu1-mu_train$mu0

  mu_test = data.frame(mu1 = bartFit$yhat.test.mean[(2*n.tr+n.ts+1):(2*n.tr+2*n.ts)],
                       mu0 = bartFit$yhat.test.mean[(2*n.tr+1):(2*n.tr+n.ts)] )
  mu_test$PLE = mu_test$mu1 - mu_test$mu0

  res = list(mods = bartFit, mu_train = mu_train, mu_test=mu_test)
  class(res) = "ple_bart"
  ## Return Results ##
  return( res )
}

#' Predict Patient-level Estimates: BART
#'
#' Get estimates of (E(Y|X,A=1), E(Y|X,A=0), E(Y|X,A=1)-E(Y|X,A=0)) using trained
#' BART model(s).
#'
#' @param object Trained BART model(s).
#' @param newdata Data-set to make predictions at.
#' @param ... Any additional parameters, not currently passed through.
#'
#' @import grf
#'
#' @return Data-frame with predictions of (E(Y|X,1), E(Y|X,0), E(Y|X,1)-E(Y|X,0))
#'
#' @examples
#' library(StratifiedMedicine)
#'
#' ## Continuous ##
#' dat_ctns = generate_subgrp_data(family="gaussian")
#' Y = dat_ctns$Y
#' X = dat_ctns$X
#' A = dat_ctns$A
#' \donttest{
#' mod1 = ple_bart(Y, A, X, Xtest=X)
#' summary(mod1$mu_train)
#' summary(predict(mod1, newdata=data.frame(A,X)))
#' summary(predict(mod1, newdata=data.frame(X)))
#' }
#' @method predict ple_bart
#' @export
#'
#### Predict: ple_bart ####
predict.ple_bart = function(object, newdata, ...){

  ### Remove Treatment variable if in newdata ##
  X = newdata[,!(colnames(newdata) %in% "A")  ]
  X.FULL = rbind( data.frame(A=0,X), data.frame(A=1,X) )
  preds = predict(object$mods, newdata=X.FULL)
  preds = apply(preds, 2, mean)
  n = dim(X)[1]
  ### PLE Predictions: Train/Test ###
  mu_hat = data.frame(mu1 = preds[(n+1):(2*n)], mu0 = preds[1:n])
  mu_hat$PLE = with(mu_hat, mu1-mu0)
  return( mu_hat  )
}
